"""The main entry point."""

import argparse
import sys
import os
from logging import getLogger


# Parse arguments FIRST so we can set the config environment variable
# BEFORE any other local modules import `config.py`.
def parse_args():
    parser = argparse.ArgumentParser(description="Run the recode_st pipeline.")
    parser.add_argument("config_file", help="Path to the TOML config file.")

    parser.add_argument(
        "--modules", nargs="+", help="List of module numbers to run (e.g. 1 2 6)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    os.environ["RECODE_CONFIG"] = args.config_file

    # Now it is safe to import local modules because config.py will read the environment variable
    from logging_config import configure_logging
    from config import settings, get_module, analysis_dir
    from seed import seed_everything
    from runtime_tracker import RuntimeTracker

    # Import pipeline steps
    from FormatData import convert_to_zarr
    from QualityControl import run_qc
    from DimensionReduction import run_dimension_reduction
    from Annotate import run_annotate
    from ViewImages import run_view_images
    from SpatialStat import run_spatial_statistics
    from MuSpan import run_muspan
    from MuSpan_SpatialGraph import run_muspan_graph
    from MuSpan_SpatialStats import run_muspan_stats
    from SelectionCSV import cosmx_csv, xenium_csv

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

    logger.info("Starting recode_st pipeline...")

    # Global settings
    data_type = settings["project"]["data_type"]
    dataset_path = settings["io"]["dataset_dir"]
    zarr_path = settings["io"]["zarr_dir"]

    # Only needed for CosMx, but won't break if missing for Xenium
    dataset_id = settings["io"].get("dataset_id", None)

    # We determine what to run based on the TOML pipeline.modules list
    if args.modules:
        valid_prefixes = {"0_", "1_", "2_", "3_", "4_", "5_", "6_"}

        modules_to_run = []
        for m in args.modules:
            prefix = f"{m}_"
            if prefix not in valid_prefixes:
                raise ValueError(f"Invalid module number: {m}")

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

            modules_to_run.append(match)
    else:
        modules_to_run = settings["pipeline"]["modules"]

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
                    qc_settings = settings["modules"]["QualityControl"]
                    run_qc(
                        data_type,
                        module_1_dir,
                        zarr_path,
                        qc_settings["min_counts"],
                        qc_settings["min_cells"],
                        qc_settings.get("min_dapi", None),  # Only for CosMx
                    )

            # MODULE 2: Dimension Reduction
            if any(m.startswith("2_") for m in modules_to_run):
                with tracker.measure("Module 2: Dimension Reduction"):
                    logger.info("Running Dimension Reduction...")
                    _, module_1_dir = get_module(1)
                    module_2_name, module_2_dir = get_module(2)
                    dr_settings = settings["modules"]["DimensionReduction"]
                    run_dimension_reduction(
                        data_type,
                        module_1_dir,
                        module_2_dir,
                        module_2_name,
                        dr_settings["n_comps"],
                        dr_settings["n_neighbors"],
                        dr_settings["resolution"],
                        dr_settings["cluster_name"],
                    )

            # MODULE 3: Annotate (Cell Type Annotation and DE Analysis)
            if any(m.startswith("3_") for m in modules_to_run):
                with tracker.measure("Module 3: Cell Annotation"):
                    logger.info("Running Annotate...")
                    _, module_2_dir = get_module(2)
                    module_3_name, module_3_dir = get_module(3)
                    cluster_name = settings["modules"]["Annotate"]["chosen_cluster"]
                    ScType_anno = settings["modules"]["Annotate"]["ScType_anno"]
                    ScType_tissue = settings["modules"]["Annotate"].get(
                        "ScType_tissue", None
                    )
                    ScType_mode = settings["modules"]["Annotate"]["ScType_mode"]
                    ScType_custom_db = settings["modules"]["Annotate"].get(
                        "ScType_custom_db", None
                    )
                    CellTypist_anno = settings["modules"]["Annotate"]["CellTypist_anno"]
                    CellTypist_model = settings["modules"]["Annotate"][
                        "CellTypist_model"
                    ]
                    CellTypist_mode = settings["modules"]["Annotate"]["CellTypist_mode"]
                    run_annotate(
                        data_type,
                        module_3_dir,
                        cluster_name,
                        module_2_dir,
                        ScType_anno,
                        ScType_tissue,
                        ScType_custom_db,
                        ScType_mode,
                        CellTypist_anno,
                        CellTypist_model,
                        CellTypist_mode,
                    )

            # MODULE 4: View Images
            if any(m.startswith("4_") for m in modules_to_run):
                with tracker.measure("Module 4: View Images"):
                    logger.info("Running View Images...")
                    viewimages_set = settings["modules"]["ViewImages"]

                    _, module_3_dir = get_module(3)
                    module_4_name, module_4_dir = get_module(4)
                    gene_list = viewimages_set["gene_list"]
                    cluster_name = settings["modules"]["Annotate"]["chosen_cluster"]

                    # xenium specific settings for grid generation in view images module
                    n_grid_x = viewimages_set.get("n_grid_x", None)
                    n_grid_y = viewimages_set.get("n_grid_y", None)

                    grid_csv_path = run_view_images(
                        data_type,
                        module_3_dir,
                        module_4_dir,
                        gene_list,
                        cluster_name,
                        n_grid_x,
                        n_grid_y,
                    )

            # MODULE 5: Spatial Statistics
            if any(m.startswith("5_") for m in modules_to_run):
                with tracker.measure("Module 5: Spatial Statistics"):
                    logger.info("Running Spatial Statistics...")
                    _, module_4_dir = get_module(4)
                    module_5_name, module_5_dir = get_module(5)
                    cluster_name = settings["modules"]["Annotate"]["chosen_cluster"]
                    run_spatial_statistics(module_5_dir, module_4_dir, cluster_name)

            # MODULE 6: MuSpAn
            if any(m.startswith("6_") for m in modules_to_run):
                with tracker.measure("Module 6: MuSpAn (Overall)"):
                    _, module_5_dir = get_module(5)
                    module_6_name, module_6_dir = get_module(6)
                    ms_settings = settings["modules"]["MuSpan"]
                    cluster_name = settings["modules"]["Annotate"]["chosen_cluster"]
                    cell_types = ms_settings["cell_types"]
                    transcript_list = ms_settings["transcripts"]

                    selection_name = ms_settings["selection_name"]
                    selected_celltypes = ms_settings["selected_celltypes"]

                    selected_fovs = ms_settings.get("selected_fovs", None)
                    box_ids = ms_settings.get("box_ids", None)

                    # --- Sub-measurements for the heaviest Module 6 parts ---
                    if data_type == "CosMx":
                        with tracker.measure("Module 6a: Create Selection CSV (CosMx)"):
                            logger.info("Creating Selection CSVs...")
                            cell_selection_csv = cosmx_csv(
                                module_dir=module_6_dir,
                                prev_module_dir=module_5_dir,
                                selection_name=selection_name,
                                cluster_col=cluster_name,
                                selected_fovs=selected_fovs,
                                selected_celltypes=selected_celltypes,
                            )

                        with tracker.measure("Module 6b: Generate MuSpAn Domain"):
                            logger.info("Running MuSpAn Domain...")
                            domain = run_muspan(
                                dataset_type=data_type,
                                module_dir=module_6_dir,
                                prev_module_dir=module_5_dir,
                                domain_name=selection_name,
                                cluster_labels=cluster_name,
                                transcripts_of_interest=transcript_list,
                                cell_selection_csv=cell_selection_csv,
                                zarr_path=zarr_path,
                                flat_files_dir=dataset_path,
                            )
                    elif data_type == "Xenium":
                        with tracker.measure(
                            "Module 6a: Create Selection CSV (Xenium)"
                        ):
                            logger.info("Creating Selection CSVs...")
                            _, module_4_dir = get_module(4)
                            grid_csv_path = (
                                f"{module_4_dir}/Xenium_ROI_grid_coordinates.csv"
                            )
                            cell_selection_csv = xenium_csv(
                                module_dir=module_6_dir,
                                prev_module_dir=module_5_dir,
                                selection_name=selection_name,
                                genes_of_interest=transcript_list,
                                cluster_col=cluster_name,
                                box_ids=box_ids,
                                grid_csv_path=grid_csv_path,
                            )

                        with tracker.measure("Module 6b: Generate MuSpAn Domain"):
                            logger.info("Running MuSpAn Domain...")
                            domain = run_muspan(
                                dataset_type=data_type,
                                module_dir=module_6_dir,
                                prev_module_dir=module_5_dir,
                                domain_name=selection_name,
                                cluster_labels=cluster_name,
                                transcripts_of_interest=transcript_list,
                                cell_selection_csv=cell_selection_csv,
                                xenium_dir=dataset_path,
                            )

                    with tracker.measure("Module 6c: MuSpAn Spatial Graphs"):
                        logger.info("Running MuSpAn Spatial Graphs...")
                        run_muspan_graph(
                            module_dir=module_6_dir,
                            domain=domain,
                            min_edge_distance=ms_settings["min_edge_distance"],
                            max_edge_distance=ms_settings["max_edge_distance"],
                            distance_list=ms_settings["distance_list"],
                            min_edge_distance_shape=ms_settings[
                                "min_edge_distance_shape"
                            ],
                            max_edge_distance_shape=ms_settings[
                                "max_edge_distance_shape"
                            ],
                            k_list=ms_settings["k_list"],
                        )

                    with tracker.measure("Module 6d: MuSpAn Spatial Stats"):
                        logger.info("Running MuSpAn Spatial Stats...")
                        run_muspan_stats(
                            module_dir=module_6_dir,
                            domain=domain,
                            cluster_labels=cluster_name,
                            cell_types=cell_types,
                        )

        logger.info(
            f"Pipeline completed successfully! Runtimes saved to: {runtime_csv_path}"
        )

    except Exception as e:
        logger.exception(f"Pipeline failed at execution: {e}")
        sys.exit(1)
