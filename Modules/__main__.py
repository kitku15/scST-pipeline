"""The main entry point."""

import argparse
import os
import sys
from logging import getLogger
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Run the recode_st pipeline.")
    parser.add_argument("config_file", help="Path to the TOML config file.")

    parser.add_argument(
        "--modules", nargs="+", help="List of module numbers to run (e.g. 1 1b 2 6)"
    )

    parser.add_argument(
        "--comp_index", type=int, help="Index of the comparison to run (0, 1, 2...)"
    )
    parser.add_argument(
        "--sample_index",
        type=int,
        help="Index of the sample to run (for PBS array jobs)",
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        help="Specific sample name to run (for manual targeting)",
    )

    return parser.parse_args()


def get_input_file(prev_dir):
    """Finds the main adata .h5ad file from the previous module."""
    prev_path = Path(prev_dir)

    # 1. Check direct path
    file_path = prev_path / "adata.h5ad"
    if file_path.exists():
        return file_path

    # 2. Check subdirectories
    sub_files = list(prev_path.glob("*/adata.h5ad"))

    if len(sub_files) == 1:
        return sub_files[0]

    elif len(sub_files) > 1:
        raise FileNotFoundError(
            f"Found multiple adata.h5ad files in {prev_dir}. "
            "You have multiple slides but skipped '1b_MergeData'. Please add it to your config!"
        )

    logger.warning(f"No adata.h5ad file found in {prev_dir}!")
    return file_path


if __name__ == "__main__":
    args = parse_args()
    os.environ["RECODE_CONFIG"] = args.config_file

    from Annotate import run_annotate
    from CellPhonedb import plot_cellphonedb, run_cellphonedb
    from config import analysis_dir, get_module, settings
    from DE_Analysis import targeted_pairwise_DE
    from Decoupler import tf_enrichment
    from DimensionReduction import run_dimension_reduction
    from FormatData import convert_to_zarr
    from LIANA_Causal import run_condition_ccc_pipeline
    from LIANA_CCC import run_liana_pipeline
    from logging_config import configure_logging
    from MergeData import run_merge
    from MuSpan import run_muspan
    from MuSpan_CellProximity import run_ms_cellproximity
    from MuSpan_Shapes import run_muspan_shapes
    from MuSpan_SpatialGraph import run_muspan_graph
    from MuSpan_SpatialStats import run_muspan_stats
    from QualityControl import run_qc
    from runtime_tracker import RuntimeTracker
    from seed import seed_everything
    from SelectionCSV import cosmx_csv, xenium_csv
    from SpatialStat import run_spatial_statistics
    from ViewImages import run_view_images
    from WebVisPrep import run_web_backend_prep

    # Setup Logging
    log_dir = analysis_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    configure_logging(log_dir)
    logger = getLogger(__package__)

    # Initialize Runtime Tracker
    runtime_csv_path = log_dir / "pipeline_runtimes.csv"
    tracker = RuntimeTracker(runtime_csv_path)

    logger.info("Seeding everything...")
    seed_everything(settings.get("seed", 42))  # Defaults to 42 if not in config

    if args.comp_index is not None:
        # Check if the module exists in the config
        if "LIANA_Causal" in settings["modules"]:
            all_comps = settings["modules"]["LIANA_Causal"].get("comparisons", [])
            if 0 <= args.comp_index < len(all_comps):
                target = all_comps[args.comp_index]
                logger.info(
                    f"PARALLEL MODE: Filtering for comparison {args.comp_index}: {target}"
                )
                # Overwrite the list with just the single targeted comparison
                settings["modules"]["LIANA_Causal"]["comparisons"] = [target]
            else:
                logger.error(
                    f"Comparison index {args.comp_index} is out of range (Total: {len(all_comps)})."
                )
                sys.exit(1)

    logger.info("Starting recode_st pipeline...")

    # Global settings
    analysis_name = settings["project"].get("analysis_name", "my_analysis")
    data_type = settings["project"]["data_type"]

    # Determine dataset path and slide name based on batch mode or single-slide mode
    base_raw_dir = settings["io"].get("base_raw_dir", None)
    
    if base_raw_dir and args.sample_index is not None:
        base_dir_path = Path(base_raw_dir)
        # Get all subdirectories, sorted alphabetically (ignores files like .tar.gz)
        dataset_folders = sorted([d for d in base_dir_path.iterdir() if d.is_dir()])
        
        # Array index is 1-based, Python lists are 0-based
        idx = args.sample_index - 1
        
        if idx >= len(dataset_folders):
            logger.info(f"Index {args.sample_index} exceeds available folders. Exiting gracefully.")
            sys.exit(0)
            
        selected_slide_dir = dataset_folders[idx]
        slide_name = selected_slide_dir.name
        
        # Override settings for downstream modules
        settings["project"]["slide_name"] = slide_name
        settings["io"]["dataset_id"] = slide_name
        dataset_path = str(selected_slide_dir)
        
        # Dynamically set the zarr path based on the folder name
        base_zarr_dir = settings["io"].get("base_zarr_dir", "data_zarrs")
        Path(base_zarr_dir).mkdir(parents=True, exist_ok=True)
        zarr_path = str(Path(base_zarr_dir) / f"{slide_name}.zarr")
        
        logger.info(f"BATCH MODE: Targeted folder '{slide_name}'.")
    else:
        # Fallback to single-slide logic if base_raw_dir isn't used
        dataset_path = settings["io"].get("dataset_dir", None)
        zarr_path = settings["io"].get("zarr_dir", None)
        slide_name = settings["project"].get("slide_name", None)

    # Batch and Sample key
    batch_key = settings["project"].get("batch_key", None)
    sample_key = settings["project"].get("sample_key", None)

    if data_type == "CosMx":
        spatial_key = "global"
    elif data_type == "Xenium":
        spatial_key = "spatial"

    # Only needed for CosMx, but won't break if missing for Xenium
    dataset_id = settings["io"].get("dataset_id", None)

    # Valid submodules for Module 6 (MuSpAn)
    valid_muspan_submodules = {"6a", "6b", "6c", "6d", "6e", "6f"}
    muspan_submodules_to_run = set()
    modules_to_run = []

    # We determine what to run based on the TOML pipeline.modules list or CLI arguments
    if args.modules:
        valid_prefixes = {
            "0_",
            "1_",
            "1b_",
            "2_",
            "3_",
            "4_",
            "5_",
            "6_",
            "7_",
            "8_",
            "8b_",
            "8c_",
            "9_",
            "10_",
        }

        for m in args.modules:
            if m == "6":
                # If "6" is selected, run all its submodules
                prefix = "6_"
                muspan_submodules_to_run.update(valid_muspan_submodules)
            elif m in valid_muspan_submodules:
                # If a specific submodule is selected (e.g. "6c")
                prefix = "6_"
                muspan_submodules_to_run.add(m)
            else:
                prefix = f"{m}_"
                if prefix not in valid_prefixes:
                    raise ValueError(f"Invalid module number or submodule: {m}")

            # Find matching full module name from config
            match = next(
                (
                    mod
                    for mod in settings["pipeline"]["modules"]
                    if mod.startswith(prefix)
                ),
                None,
            )

            if match is None:
                raise ValueError(f"No module found in config for prefix {prefix}")

            # Prevent duplicate appending if multiple 6 submodules are called
            if match not in modules_to_run:
                modules_to_run.append(match)
    else:
        modules_to_run = settings["pipeline"]["modules"]
        # If running purely from config and Module 6 is included, run all its submodules by default
        if any(m.startswith("6_") for m in modules_to_run):
            muspan_submodules_to_run.update(valid_muspan_submodules)

    try:
        with tracker.measure("Total Pipeline Execution"):
            # MODULE 0: Format Data
            if any(m.startswith("0_") for m in modules_to_run):
                with tracker.measure("Module 0: Format Data"):
                    logger.info("Running Format Data...")
                    convert_to_zarr(data_type, dataset_path, dataset_id, zarr_path)

            # MODULE 1: Quality Control
            if any(m.startswith("1_") for m in modules_to_run):
                with tracker.measure("Module 1: Quality Control"):
                    logger.info("Running Quality Control...")
                    module_1_name, module_1_dir = get_module(1)

                    slide_name = settings["project"].get("slide_name", None)
                    if slide_name:
                        module_1_dir = module_1_dir / slide_name
                        module_1_dir.mkdir(parents=True, exist_ok=True)

                    qc_settings = settings["modules"]["QualityControl"]

                    # replace {slide_name} placeholder
                    fov_path = qc_settings.get("fov_metadata_path", None)
                    if fov_path and "{slide_name}" in fov_path:
                        fov_path = fov_path.format(slide_name=slide_name)
                        
                    proseg_path = qc_settings.get("proseg_zarr_path", None)
                    if proseg_path and "{slide_name}" in proseg_path:
                        proseg_path = proseg_path.format(slide_name=slide_name)

                    run_qc(
                        data_type=data_type,
                        module_dir=module_1_dir,
                        zarr_path=zarr_path,
                        min_counts=qc_settings["min_counts"],
                        min_cells=qc_settings["min_cells"],
                        min_genes=qc_settings["min_genes"],
                        min_area=qc_settings["min_area"],
                        max_area=qc_settings["max_area"],
                        min_dapi=qc_settings.get("min_dapi", None),
                        batch_key=batch_key,
                        sample_key=sample_key,
                        fov_metadata_path=fov_path,
                        proseg_zarr_path=proseg_path,
                        proseg_cell_id_col=qc_settings.get(
                            "proseg_cell_id_col", "original_cell_id"
                        ),
                        original_cell_id_col=qc_settings.get(
                            "original_cell_id_col", "index"
                        ),
                    )

            # MODULE 1B: Merge Data
            if any(m.startswith("1b_") for m in modules_to_run):
                with tracker.measure("Module 1b: Merge Data"):
                    logger.info("Running Data Merging...")
                    module_1b_name, module_1b_dir = get_module("1b")
                    _, module_1_dir = get_module(1)
                    
                    merge_settings = settings["modules"].get("MergeData", {})

                    # Try to get them from config, otherwise AUTO-DISCOVER them
                    input_files = merge_settings.get("input_files", [])
                    slide_names = merge_settings.get("slide_names", [])

                    if not input_files:
                        logger.info(f"Auto-discovering QC'd datasets in {module_1_dir}...")
                        # Find all adata.h5ad files in the subdirectories of Module 1
                        found_files = sorted(list(module_1_dir.glob("*/adata.h5ad")))
                        
                        if not found_files:
                            raise FileNotFoundError(f"No QC'd adata.h5ad files found in {module_1_dir}. Did Module 1 finish?")
                            
                        input_files = [str(f) for f in found_files]
                        # The slide name is the name of the folder containing the adata.h5ad
                        slide_names = [f.parent.name for f in found_files]
                        
                        logger.info(f"Found {len(input_files)} datasets to merge: {slide_names}")

                    run_merge(
                        input_files=input_files,
                        slide_names=slide_names,
                        module_dir=module_1b_dir,
                    )

            # MODULE 2: Dimension Reduction
            if any(m.startswith("2_") for m in modules_to_run):
                with tracker.measure("Module 2: Dimension Reduction"):
                    logger.info("Running Dimension Reduction...")

                    # Dynamically look for input: If 1b was defined, use it. Otherwise fallback to 1.
                    try:
                        _, prev_dir = get_module("1b")
                    except ValueError:
                        _, prev_dir = get_module(1)

                    module_2_name, module_2_dir = get_module(2)
                    dr_settings = settings["modules"]["DimensionReduction"]

                    input_file = get_input_file(prev_dir)
                    run_dimension_reduction(
                        data_type=data_type,
                        input_adata_path=input_file,
                        module_dir=module_2_dir,
                        module_name=module_2_name,
                        n_comps=dr_settings.get("n_comps", None),
                        use_scviva=dr_settings.get("use_scviva", True),
                        n_neighbors_list=dr_settings["n_neighbors"],
                        resolution_list=dr_settings["resolution"],
                        cluster_name=dr_settings["cluster_name"],
                        umap_latent=dr_settings["umap_latent"],
                        scviva_layer=dr_settings["scviva_layer"],
                        scviva_batch_key=batch_key,
                        scviva_sample_key=sample_key,
                        scviva_spatial_knn=dr_settings["scviva_spatial_knn"],
                        scviva_epochs=dr_settings["scviva_epochs"],
                        run_pca=dr_settings.get("run_pca", False),
                        use_scanvi=dr_settings.get("use_scanvi", False),
                        scanvi_mode=dr_settings.get("scanvi_mode", "scarches"),
                        reference_adata_path=dr_settings.get(
                            "reference_adata_path", None
                        ),
                        reference_batch_key=dr_settings.get(
                            "reference_batch_key", None
                        ),
                        reference_label_key=dr_settings.get(
                            "reference_label_key", None
                        ),
                        reference_layer=dr_settings.get("reference_layer", "X"),
                        scvi_epochs=dr_settings.get("scvi_epochs", 400),
                        scanvi_epochs=dr_settings.get("scanvi_epochs", 200),
                        scviva_batch_size=dr_settings.get("scviva_batch_size", 512),
                        n_latent=dr_settings.get("n_latent", 30),
                        n_layers=dr_settings.get("n_layers", 2),
                        pre_cluster_res=dr_settings.get("pre_cluster_res", 1.0),
                        dot_size=dr_settings.get("dot_size", 5),
                    )

            # MODULE 3: Annotate (Cell Type Annotation and DE Analysis)
            if any(m.startswith("3_") for m in modules_to_run):
                with tracker.measure("Module 3: Cell Annotation"):
                    logger.info("Running Annotate...")
                    _, module_2_dir = get_module(2)
                    module_3_name, module_3_dir = get_module(3)
                    anno_settings = settings["modules"]["Annotate"]

                    input_file = get_input_file(module_2_dir)
                    run_annotate(
                        datatype=data_type,
                        module_dir=module_3_dir,
                        cluster_name=anno_settings["chosen_cluster"],
                        input_adata_path=input_file,
                        sample_key=sample_key,
                        ScType_anno=anno_settings.get("ScType_anno", False),
                        ScType_tissue=anno_settings.get("ScType_tissue", None),
                        ScType_custom_db=anno_settings.get("ScType_custom_db", None),
                        ScType_mode=anno_settings.get("ScType_mode", None),
                        CellTypist_anno=anno_settings["CellTypist_anno"],
                        CellTypist_model=anno_settings.get("CellTypist_model", None),
                        CellTypist_mode=anno_settings["CellTypist_mode"],
                        CellTypist_custom_model=anno_settings.get(
                            "CellTypist_custom_model", None
                        ),
                        CellTypist_train=anno_settings.get("CellTypist_train", False),
                        CellTypist_train_data=anno_settings.get(
                            "CellTypist_train_data", None
                        ),
                        CellTypist_train_labels=anno_settings.get(
                            "CellTypist_train_labels", None
                        ),
                        plot=anno_settings.get("plot", True),
                    )

            # MODULE 4: View Images
            if any(m.startswith("4_") for m in modules_to_run):
                with tracker.measure("Module 4: View Images"):
                    logger.info("Running View Images...")
                    _, module_3_dir = get_module(3)
                    module_4_name, module_4_dir = get_module(4)
                    viewimages_set = settings["modules"]["ViewImages"]

                    input_file = get_input_file(module_3_dir)
                    grid_csv_path = run_view_images(
                        data_type=data_type,
                        input_adata_path=input_file,
                        sample_key=sample_key,
                        module_dir=module_4_dir,
                        gene_list=viewimages_set["gene_list"],
                        cluster_name=viewimages_set["chosen_cluster"],
                        n_grid_x=viewimages_set.get("n_grid_x", 10),
                        n_grid_y=viewimages_set.get("n_grid_y", 10),
                        embedding_key=viewimages_set.get("embedding_key", None),
                        umap_color_columns=viewimages_set.get(
                            "umap_color_columns", None
                        ),
                    )

            # MODULE 5: Spatial Statistics
            if any(m.startswith("5_") for m in modules_to_run):
                with tracker.measure("Module 5: Spatial Statistics"):
                    logger.info("Running Spatial Statistics...")
                    _, module_4_dir = get_module(4)
                    module_5_name, module_5_dir = get_module(5)

                    spatial_settings = settings["modules"]["SpatialStat"]
                    cluster_name = spatial_settings["chosen_cluster"]
                    condition_key = spatial_settings.get("condition_key", None)
                    reference_condition = spatial_settings.get(
                        "reference_condition", None
                    )

                    input_file = get_input_file(module_4_dir)
                    run_spatial_statistics(
                        module_dir=module_5_dir,
                        input_adata_path=input_file,
                        sample_key=sample_key,
                        cluster_name=cluster_name,
                        condition_key=condition_key,
                        reference_condition=reference_condition,
                        skip_compute=False,
                    )

            # MODULE 6: MuSpAn
            if any(m.startswith("6_") for m in modules_to_run):
                with tracker.measure("Module 6: MuSpAn (Overall)"):
                    _, module_5_dir = get_module(5)
                    module_6_name, module_6_dir = get_module(6)
                    ms_settings = settings["modules"]["MuSpan"]
                    cell_types = ms_settings["cell_types"]
                    transcript_list = ms_settings["transcripts"]
                    base_selection_name = ms_settings.get("selection_name", None)
                    cluster_name = ms_settings["chosen_cluster"]
                    selected_celltypes = ms_settings.get("selected_celltypes", None)

                    # FETCH SELECTED FOVS AS A DICTIONARY
                    selected_fovs_dict = ms_settings.get("selected_fovs", {})

                    input_file = get_input_file(module_5_dir)

                    import scanpy as sc

                    logger.info("Loading merged adata to process MuSpAn per sample...")
                    adata_merged = sc.read_h5ad(input_file)

                    # FILTER SAMPLES BASED ON CLI ARGUMENTS
                    all_samples = sorted(list(adata_merged.obs[sample_key].unique()))
                    samples_to_run = all_samples

                    if args.sample_index is not None:
                        samples_to_run = [all_samples[args.sample_index]]
                    elif args.sample_name is not None:
                        samples_to_run = [args.sample_name]

                    # GET LIST OF SELECTIONS FROM TOML
                    base_selections = ms_settings.get("selection_names", None)

                    # 3. loop throug the Samples
                    for sample in samples_to_run:
                        sample_out_dir = module_6_dir / sample
                        sample_out_dir.mkdir(parents=True, exist_ok=True)

                        sample_input_file = (
                            sample_out_dir / f"adata_slice_{sample}.h5ad"
                        )
                        if not sample_input_file.exists():
                            adata_slice = adata_merged[
                                adata_merged.obs[sample_key] == sample
                            ].copy()
                            adata_slice.write_h5ad(sample_input_file)

                        ## 4. loop through the Regions (e.g., myeloid_ulcer, healthy_ulcer)
                        for base_selection_name in base_selections:
                            logger.info(
                                f"=== Processing MuSpAn for Sample: {sample} | Region: {base_selection_name} ==="
                            )

                            selection_name = f"{base_selection_name}_{sample}"
                            selection_out_dir = sample_out_dir / selection_name
                            selection_out_dir.mkdir(parents=True, exist_ok=True)

                            # Only skip if actual MuSpAn outputs (like images) exist
                            is_done = False
                            if selection_out_dir.exists():
                                # Check if there are files OTHER than just the CSV
                                files = [
                                    f
                                    for f in selection_out_dir.iterdir()
                                    if not f.name.endswith(".csv")
                                ]
                                if len(files) > 0:
                                    is_done = True

                            if is_done:
                                logger.info(
                                    f"Outputs for {selection_name} already exist. Skipping..."
                                )
                                continue

                            slide_id = adata_merged.obs[
                                adata_merged.obs[sample_key] == sample
                            ]["slide_id"].iloc[0]

                            raw_io = (
                                settings["io"].get("raw_data", {}).get(slide_id, {})
                            )
                            dataset_path = raw_io.get("dataset_dir")
                            zarr_path = raw_io.get("zarr_dir")
                            proseg_zarr_path = raw_io.get("proseg_zarr_dir")

                            if not dataset_path:
                                logger.error(
                                    f"Missing [io.raw_data] config for {slide_id}. Skipping MuSpAn for {sample}."
                                )
                                continue

                            sample_fovs = selected_fovs_dict.get(sample, None)

                            # get colors for muspan plotting
                            color_dict = None
                            colors_key = f"{cluster_name}_colors"
                            if (
                                colors_key in adata_merged.uns
                                and cluster_name in adata_merged.obs.columns
                            ):
                                categories = list(
                                    adata_merged.obs[cluster_name].cat.categories
                                )
                                colors = list(adata_merged.uns[colors_key])
                                color_dict = dict(zip(categories, colors))

                            domain = None
                            cell_selection_csv = None

                            # Check for Custom Pre-Generated CSV
                            expected_csv = selection_out_dir / f"{selection_name}.csv"

                            if expected_csv.exists():
                                cell_selection_csv = expected_csv

                            if cell_selection_csv:
                                logger.info(
                                    f"Custom CSV found! Bypassing Module 6a and using: {cell_selection_csv}"
                                )

                            # MODULE 6A: Create Selection CSV (Skipped if custom CSV is found)
                            if (
                                "6a" in muspan_submodules_to_run
                                and cell_selection_csv is None
                            ):
                                if data_type == "CosMx":
                                    with tracker.measure(
                                        f"Module 6a: Selection CSV (CosMx) - {selection_name}"
                                    ):
                                        cell_selection_csv = cosmx_csv(
                                            module_dir=sample_out_dir,
                                            input_adata_path=sample_input_file,
                                            selection_name=selection_name,
                                            cluster_col=cluster_name,
                                            selected_fovs=sample_fovs,
                                            selected_celltypes=selected_celltypes,
                                        )
                                elif data_type == "Xenium":
                                    with tracker.measure(
                                        f"Module 6a: Selection CSV (Xenium) - {selection_name}"
                                    ):
                                        _, module_4_dir = get_module(4)
                                        grid_csv_path = (
                                            module_4_dir
                                            / f"Xenium_ROI_grid_coordinates_{sample}.csv"
                                        )
                                        cell_selection_csv = xenium_csv(
                                            module_dir=sample_out_dir,
                                            input_adata_path=sample_input_file,
                                            selection_name=selection_name,
                                            genes_of_interest=transcript_list,
                                            cluster_col=cluster_name,
                                            box_ids=sample_fovs,
                                            grid_csv_path=grid_csv_path,
                                        )

                            # MODULE 6B: Generate MuSpAn Domain
                            if "6b" in muspan_submodules_to_run:
                                if data_type == "CosMx":
                                    with tracker.measure(
                                        f"Module 6b: Generate MuSpAn Domain - {selection_name}"
                                    ):
                                        domain = run_muspan(
                                            dataset_type=data_type,
                                            module_dir=sample_out_dir,
                                            input_adata_path=sample_input_file,
                                            domain_name=selection_name,
                                            cluster_labels=cluster_name,
                                            transcripts_of_interest=transcript_list,
                                            cell_selection_csv=cell_selection_csv,
                                            zarr_path=zarr_path,
                                            flat_files_dir=dataset_path,
                                            proseg_zarr_path=proseg_zarr_path,
                                            color_dict=color_dict,
                                        )
                                elif data_type == "Xenium":
                                    with tracker.measure(
                                        f"Module 6b: Generate MuSpAn Domain - {selection_name}"
                                    ):
                                        domain = run_muspan(
                                            dataset_type=data_type,
                                            module_dir=sample_out_dir,
                                            input_adata_path=sample_input_file,
                                            domain_name=selection_name,
                                            cluster_labels=cluster_name,
                                            transcripts_of_interest=transcript_list,
                                            cell_selection_csv=cell_selection_csv,
                                            xenium_dir=dataset_path,
                                        )

                            # MODULE 6C: MuSpAn Spatial Graphs
                            if "6c" in muspan_submodules_to_run:
                                with tracker.measure(
                                    f"Module 6c: MuSpAn Spatial Graphs - {selection_name}"
                                ):
                                    run_muspan_graph(
                                        module_dir=sample_out_dir,
                                        domain=domain,
                                        min_edge_distance=ms_settings[
                                            "min_edge_distance"
                                        ],
                                        max_edge_distance=ms_settings[
                                            "max_edge_distance"
                                        ],
                                        distance_list=ms_settings["distance_list"],
                                        min_edge_distance_shape=ms_settings[
                                            "min_edge_distance_shape"
                                        ],
                                        max_edge_distance_shape=ms_settings[
                                            "max_edge_distance_shape"
                                        ],
                                        k_list=ms_settings["k_list"],
                                    )

                            # MODULE 6D: MuSpAn Spatial Statistics
                            if "6d" in muspan_submodules_to_run:
                                with tracker.measure(
                                    f"Module 6d: MuSpAn Spatial Stats - {selection_name}"
                                ):
                                    run_muspan_stats(
                                        module_dir=sample_out_dir,
                                        domain=domain,
                                        cluster_labels=cluster_name,
                                        cell_types=cell_types,
                                    )

                            # MODULE 6E: MuSpAn Shape Analysis
                            if "6e" in muspan_submodules_to_run:
                                with tracker.measure(
                                    f"Module 6e: MuSpAn Shape Analysis - {selection_name}"
                                ):
                                    run_muspan_shapes(
                                        module_dir=sample_out_dir,
                                        domain=domain,
                                        chosen_cluster=cluster_name,
                                        selected_celltypes=selected_celltypes,
                                    )

                            # MODULE 6F: MuSpAn Cell Proximity Analysis
                            if "6f" in muspan_submodules_to_run:
                                with tracker.measure(
                                    f"Module 6f: MuSpAn Cell Proximity Analysis - {selection_name}"
                                ):
                                    run_ms_cellproximity(
                                        module_dir=sample_out_dir,
                                        domain=domain,
                                        chosen_cluster=cluster_name,
                                        selected_celltypes=ms_settings[
                                            "network_celltypes"
                                        ],
                                        selection_name=selection_name,
                                        max_distance=ms_settings.get(
                                            "max_distance", 20
                                        ),
                                        cellboundary_label=ms_settings.get(
                                            "cellboundary_label", "Cell boundaries"
                                        ),
                                        color_dict=color_dict,
                                    )
                    # MUSPAN CONDITION AGGREGATION
                    # Only run if condition_key is defined in config
                    condition_key = ms_settings.get("condition_key", None)
                    reference_condition = ms_settings.get("reference_condition", None)

                    if condition_key and condition_key in adata_merged.obs.columns:
                        if args.sample_index is None and args.sample_name is None:
                            sample_to_cond = (
                                adata_merged.obs.groupby(sample_key)[condition_key]
                                .first()
                                .to_dict()
                            )

                            from MuSpan_SpatialStats import (
                                aggregate_muspan_across_conditions,
                            )

                            with tracker.measure(
                                "Module 6: Cross-Condition Aggregation"
                            ):
                                aggregate_muspan_across_conditions(
                                    module_6_dir=module_6_dir,
                                    sample_condition_dict=sample_to_cond,
                                    reference_condition=reference_condition,
                                    base_selection_name=base_selections[0]
                                    if "base_selections" in locals()
                                    else ms_settings.get("selection_name"),
                                )
                        else:
                            logger.info(
                                "Skipping condition aggregation because a specific sample was targeted."
                            )
                    else:
                        logger.warning(
                            "No 'condition_key' found in [modules.MuSpan]. Skipping MuSpAn condition aggregation."
                        )

            if any(m.startswith("7_") for m in modules_to_run):
                with tracker.measure("Module 7: Decoupler TF Enrichment Analysis"):
                    logger.info("Running Decoupler TF Enrichment Analysis...")
                    _, module_7_dir = get_module(7)
                    _, module_5_dir = get_module(5)
                    _, module_10_dir = get_module(10)
                    decoupler_settings = settings["modules"]["Decoupler"]
                    organism = decoupler_settings.get("organism", None)
                    grn = decoupler_settings.get("grn", None)
                    dorothea_levels = decoupler_settings.get("dorothea_levels", None)
                    celltype_key = decoupler_settings.get("celltype_key", None)
                    active_tfs_file_name = decoupler_settings.get(
                        "active_tfs_file_name", "active_tf.txt"
                    )

                    input_file = get_input_file(module_5_dir)
                    tf_enrichment(
                        module_dir=module_7_dir,
                        web_dir=module_10_dir,
                        input_adata_path=input_file,
                        sample_key=sample_key,
                        celltype_key=celltype_key,
                        organism=organism,
                        grn=grn,
                        dorothea_levels=["A", "B", "C"],
                        active_tfs_file_name=active_tfs_file_name,
                    )

            if any(m.startswith("8_") for m in modules_to_run):
                _, module_8_dir = get_module(8)
                _, module_5_dir = get_module(5)
                _, module_7_dir = get_module(7)
                cpdb_settings = settings["modules"]["Cellphonedb"]
                chosen_cluster = cpdb_settings.get("chosen_cluster", None)
                cpdb_version = cpdb_settings.get("cpdb_version", "v5.0.0")
                human = cpdb_settings.get("human", True)
                celltypes = cpdb_settings.get("celltypes", None)
                target_genes = cpdb_settings.get("target_genes", None)
                gene_family = cpdb_settings.get("gene_family", None)
                microenv_res = cpdb_settings.get("microenv_resolution", 0.5)
                target_microenv = cpdb_settings.get("target_microenv", None)
                target_sample = cpdb_settings.get("target_sample", None)
                use_cellsign = cpdb_settings.get("use_cellsign", False)

                decoupler_settings = settings["modules"]["Decoupler"]
                active_tfs_file_name = decoupler_settings.get(
                    "active_tfs_file_name", "active_tf.txt"
                )
                active_tfs_file_path = module_7_dir / active_tfs_file_name

                with tracker.measure("Module 8: Cellphonedb CCC Analysis"):
                    logger.info("Running Cellphonedb CCC Analysis...")
                    input_file = get_input_file(module_5_dir)

                    cpdb_counts_path = run_cellphonedb(
                        module_dir=module_8_dir,
                        input_adata_path=input_file,
                        sample_key=sample_key,
                        chosen_cluster=chosen_cluster,
                        cpdb_version=cpdb_version,
                        active_tfs_file_path=active_tfs_file_path,
                        human=human,
                        microenv_resolution=microenv_res,
                    )
                    logger.info("Plotting Cellphonedb Results...")

                    plot_cellphonedb(
                        module_dir=module_8_dir,
                        cpdb_counts_path=cpdb_counts_path,
                        celltype_key=chosen_cluster,
                        celltypes=celltypes,
                        target_microenv=target_microenv,
                        target_sample=target_sample,
                        sample_key=sample_key,
                        gene_family=gene_family,
                        target_genes=target_genes,
                        use_cellsign=use_cellsign,
                    )

            # MODULE 8b: LIANA+ Spatial CCC
            if any(m.startswith("8b_") for m in modules_to_run):
                _, module_8b_dir = get_module("8b")
                _, module_5_dir = get_module(5)
                _, module_7_dir = get_module(7)

                liana_settings = settings["modules"].get("LIANA", {})
                tf_file = module_7_dir / "tf_activity_scores.h5ad"

                with tracker.measure("Module 8b: LIANA+ Single-Cell CCC"):
                    logger.info("Running LIANA+ Spatial CCC...")
                    input_file = get_input_file(module_5_dir)

                    run_liana_pipeline(
                        module_dir=module_8b_dir,
                        input_adata_path=input_file,
                        tf_adata_path=tf_file,
                        spatial_key=spatial_key,
                        sample_key=sample_key,
                        settings=liana_settings,
                    )

            # MODULE 8c: LIANA+ Condition-Specific CCC & Causal Inference
            if any(m.startswith("8c_") for m in modules_to_run):
                _, module_8c_dir = get_module("8c")
                _, module_5_dir = get_module(5)
                _, module_7_dir = get_module(7)

                liana_causal_settings = settings["modules"].get("LIANA_Causal", {})

                with tracker.measure(
                    "Module 8c: LIANA+ Condition CCC & Causal Inference"
                ):
                    logger.info("Running LIANA+ Condition-Specific Causal Inference...")
                    input_file = get_input_file(module_5_dir)
                    pseudobulk_file = module_7_dir / "pseudobulk.h5ad"

                    run_condition_ccc_pipeline(
                        module_dir=module_8c_dir,
                        input_adata_path=input_file,
                        pseudobulk_adata_path=pseudobulk_file,
                        settings=liana_causal_settings,
                    )

            if any(m.startswith("9_") for m in modules_to_run):
                _, module_9_dir = get_module(9)
                _, module_7_dir = get_module(7)
                DEAnalysis_settings = settings["modules"]["DEAnalysis"]

                with tracker.measure("Module 9: Pydeseq2 Pseudobulk Analysis"):
                    targeted_pairwise_DE(
                        pseudobulk_adata_path=f"{module_7_dir}/pseudobulk.h5ad",
                        celltype_col=DEAnalysis_settings.get("celltype_col", None),
                        treatment_col=DEAnalysis_settings.get("treatment_col", None),
                        comparisons=DEAnalysis_settings.get("comparisons", None),
                        module_dir=module_9_dir,
                    )

            # MODULE 10: Exporting Files for Web Tool
            if any(m.startswith("10_") for m in modules_to_run):
                _, module_1_dir = get_module(1)
                _, module_3_dir = get_module(3)
                _, module_5_dir = get_module(5)
                _, module_6_dir = get_module(6)
                _, module_7_dir = get_module(7)
                _, module_8_dir = get_module(8)
                _, module_8b_dir = get_module("8b")
                _, module_8c_dir = get_module("8c")
                _, module_9_dir = get_module(9)
                _, module_10_dir = get_module(10)
                
                # Fetch the heavily processed pseudobulk/TF AnnData
                input_file = get_input_file(module_7_dir)

                DEAnalysis = settings["modules"]["DEAnalysis"].get("DEAnalysis", False)
                webvissettings = settings["modules"]["WebVisPrep"]

                celltype_key = webvissettings.get("primary_annotation", "Broad_Celltype")
                microenv_key = webvissettings.get("microenv_col", "spatial_microenvironment")
                anno_keywords = webvissettings.get("annotation_columns", ["leiden", "CellTypist", "sctype", "cluster"])

                with tracker.measure("Module 10: Exporting Zarr and Aux Data for Web Backend"):
                    
                    run_web_backend_prep(
                        module_dir=module_10_dir,
                        input_adata_path=input_file,
                        batch_key=batch_key,
                        sample_key=sample_key,
                        celltype_key=celltype_key,
                        microenv_key=microenv_key,
                        module_1_dir=module_1_dir,
                        module_3_dir=module_3_dir,
                        module_5_dir=module_5_dir,
                        module_6_dir=module_6_dir,
                        module_7_dir=module_7_dir,
                        module_8_dir=module_8_dir,
                        module_8b_dir=module_8b_dir,
                        module_8c_dir=module_8c_dir,
                        module_9_dir=module_9_dir,
                        analysis_name=analysis_name,
                        spatial_key=spatial_key,
                        data_type=data_type,
                        settings=settings,
                        DEAnalysis=DEAnalysis,
                        anno_keywords=anno_keywords
                    )

        logger.info(
            f"Pipeline completed successfully! Runtimes saved to: {runtime_csv_path}"
        )

    except Exception as e:
        logger.exception(f"Pipeline failed at execution: {e}")
        sys.exit(1)
