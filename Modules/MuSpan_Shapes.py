"""Muspan module - Shapes analysis."""

import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
logger = getLogger(__name__)

try:
    import muspan as ms
except ModuleNotFoundError as err:
    logger.info(
        "Could not load MuSpAn. Install with:\n"
        "    pip install 'recode_st[muspan]' @ git+https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
    )
    raise err


def ms_shapedesc(domain):
    # Compute metrics
    logger.info("Calculating shape descriptors...")

    # 1. Area
    ms.geometry.area(
        domain,
        population=("collection", "Cell boundaries"),
        add_as_label=True,
        label_name="Area (µm²)",
        cmap="inferno",
    )
    # 2. Perimeter
    ms.geometry.perimeter(
        domain,
        population=("collection", "Cell boundaries"),
        add_as_label=True,
        label_name="Perimeter (µm)",
        cmap="cividis",
    )
    # 3. Convexity
    ms.geometry.convexity(
        domain,
        population=("collection", "Cell boundaries"),
        add_as_label=True,
        label_name="Convexity",
        cmap="bone",
    )
    # 4. Circularity
    ms.geometry.circularity(
        domain,
        population=("collection", "Cell boundaries"),
        add_as_label=True,
        label_name="Circularity",
        cmap="viridis",
    )

    # Create 2x2 subplot figure
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    plot_info = [
        ("Area (µm²)", "inferno", axes[0, 0]),
        ("Perimeter (µm)", "cividis", axes[0, 1]),
        ("Convexity", "bone", axes[1, 0]),
        ("Circularity", "viridis", axes[1, 1]),
    ]

    for label, cmap, ax in plot_info:
        plt.sca(ax)

        ms.visualise.visualise(
            domain,
            color_by=label,
            objects_to_plot=("collection", "Cell boundaries"),
            shape_kwargs=dict(alpha=1),
            ax=ax,
        )

        ax.set_title(label)

    return plt


def ms_PrAxis(domain, chosen_cluster, cluster_list=None):
    logger.info("Calculating shape orientation via principle axis...")

    if cluster_list is None:
        # Retrieve all values for the chosen_cluster label
        labels_array, object_indices = ms.query.get_labels(domain, chosen_cluster)
        # Extract the unique values
        cluster_list = np.unique(labels_array)

    target_cells = ms.query.query_container(
        ("collection", "Cell boundaries"),
        "AND",
        ms.query.query(domain, ("label", chosen_cluster), "in", cluster_list),
        domain,
    )

    # 2. Calculate angles
    angles_dict = {}
    for cluster in cluster_list:
        cluster_cells = ms.query.query_container(
            ("collection", "Cell boundaries"), "AND", (chosen_cluster, cluster), domain
        )
        angles, _, _ = ms.geometry.principle_axis(
            domain,
            population=cluster_cells,
            add_as_label=True,
            label_name="Principle axis angle (rad)",
            cmap="coolwarm",
        )
        angles_dict[cluster] = angles

    # 3. Set up the figure
    fig = plt.figure(figsize=(10, 7))
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    for cluster in sorted([str(c) for c in cluster_list]):
        angles = angles_dict[str(cluster)]  # Pull the correct angles
        sns.kdeplot(
            angles,
            ax=ax1,
            label=f"{cluster} Cells",
            fill=True,
            clip=(-3.14 / 2, 3.14 / 2),
        )

    ax1.set_xlabel("Principle axis angle (rad)")
    ax1.set_ylabel("Density")
    ax1.legend()

    # 5. Visualise domain with principle axis
    ms.visualise.visualise(
        domain,
        color_by=("constant", [0.7, 0.7, 0.7, 1]),
        objects_to_plot=("collection", "Cell boundaries"),
        shape_kwargs=dict(alpha=0.5, linewidth=0.5),
        ax=ax3,
    )
    ms.visualise.visualise(
        domain,
        color_by="Principle axis angle (rad)",
        objects_to_plot=target_cells,
        shape_kwargs=dict(alpha=1, linewidth=0.5),
        ax=ax3,
        vmin=-3.14 / 2,
        vmax=3.14 / 2,
    )

    # 6. EXACT ORIGINAL CODE: Visualise the chosen_cluster
    # Because we aren't interfering, MuSpAn will automatically draw its perfect native colorbar
    ms.visualise.visualise(
        domain,
        color_by=("constant", [0.7, 0.7, 0.7, 1]),
        objects_to_plot=("collection", "Cell boundaries"),
        shape_kwargs=dict(alpha=0.5, linewidth=0.5),
        ax=ax2,
    )
    ms.visualise.visualise(
        domain,
        color_by=chosen_cluster,
        objects_to_plot=target_cells,
        shape_kwargs=dict(alpha=1, linewidth=0.5),
        ax=ax2,
    )

    return plt


def ms_PointstoShape(domain, chosen_cluster, cell_type):
    # 1. Query the specific Leiden cluster (replace '1' with your target cluster)
    cluster_cells = ms.query.query(domain, ("label", chosen_cluster), "is", cell_type)

    # 2. Query to ensure we are only looking at Cell centroids
    centroids = ms.query.query(domain, ("collection",), "is", "Cell centroids")

    # 3. Combine them so you only get Cell centroids that belong to cluster '1'
    target_population = ms.query.query_container(cluster_cells, "AND", centroids)

    # Visualize to ensure you've grabbed the right points before proceeding
    # ms.visualise.visualise(
    #     domain, color_by=("label", chosen_cluster), objects_to_plot=target_population
    # )

    # Create a figure with 4 subplots
    fig, axes = plt.subplots(1, 4, figsize=(12, 3))

    # 1. Try significantly larger alpha values based on spatial transcriptomics scales
    # You may need to add a zero to these if your dataset is very large!
    test_alphas = [30, 50, 100, 200]

    for i, alpha in enumerate(test_alphas):
        try:
            # Convert the objects to a single shape using the alpha shape method
            new_IDs = domain.convert_objects(
                population=target_population,
                object_type="shape",
                conversion_method="alpha shape",
                conversion_method_kwargs=dict(alpha=alpha),
                collection_name=f"Alpha Object {alpha}",
                inherit_collections=False,
                return_IDs=True,
            )

            # Visualize the original cells in grey
            ms.visualise.visualise(
                domain,
                color_by=("constant", "grey"),
                objects_to_plot=target_population,
                ax=axes[i],
                marker_size=1,
            )

            # Visualize the new alpha shape in red
            ms.visualise.visualise(
                domain,
                color_by=("constant", "red"),
                objects_to_plot=new_IDs,
                ax=axes[i],
            )

            # Set the title for successful subplots
            axes[i].set_title(f"Alpha = {alpha}")

        except ValueError:
            # 2. If the alpha is too small and generates the zero-size array error, catch it
            print(f"Skipping Alpha = {alpha}: Value too small to connect any points.")

            # Plot just the points so you can see why it failed
            ms.visualise.visualise(
                domain,
                color_by=("constant", "grey"),
                objects_to_plot=target_population,
                ax=axes[i],
                marker_size=1,
            )
            axes[i].set_title(f"Alpha = {alpha}\n(Failed - Too Small)")

    plot_name = "ms_pts.png"
    plt.tight_layout()
    plt.savefig(plot_name)
    logger.info(f"points to shape plot successfully saved as {plot_name}")


def run_muspan_shapes(module_dir, domain, chosen_cluster, selected_celltypes):
    domain_name = domain.name

    # 1. Shape descriptors
    plt = ms_shapedesc(domain)

    out_dir = Path(module_dir) / domain_name
    plot_save = f"{out_dir}/ms_shapes.png"
    plt.tight_layout()
    plt.savefig(plot_save)
    logger.info(f"Shape descriptor plot successfully saved as {plot_save}")
    plt.close()

    # 2. Principle axis
    plt = ms_PrAxis(domain, chosen_cluster, selected_celltypes)
    out_dir = Path(module_dir) / domain_name
    plot_save = f"{out_dir}/ms_praxis.png"
    plt.tight_layout()
    plt.savefig(plot_save)
    logger.info(f"shape orientation plot successfully saved as {plot_save}")

    # Export Morphometrics to CSV
    logger.info("Extracting shape metrics into a CSV for web visualization...")

    try:
        # 1. Grab the "Cell ID" labels directly (Bypasses collection queries!)
        # This returns numpy arrays which are safely iterable
        cell_ids, cell_id_obj_indices = ms.query.get_labels(domain, "Cell ID")
        cell_id_dict = dict(zip(list(cell_id_obj_indices), list(cell_ids)))

        # 2. Grab the "Area (µm²)" labels to identify our boundaries
        # Because Area is only calculated on boundaries, this perfectly isolates them
        areas, boundary_indices = ms.query.get_labels(domain, "Area (µm²)")
        boundary_list = list(boundary_indices)

        logger.info(
            f"DEBUG: Found {len(cell_id_dict)} Cell IDs and {len(boundary_list)} Boundaries."
        )

        # 3. Create DataFrame
        df_morph = pd.DataFrame({"boundary_id": boundary_list})

        # 4. Safely map Cell IDs to Boundaries
        # CASE A: Xenium (The boundaries themselves hold the Cell ID string)
        if len(boundary_list) > 0 and boundary_list[0] in cell_id_dict:
            logger.info("DEBUG: Mapping Cell IDs directly from boundaries.")
            df_morph["Cell_ID"] = df_morph["boundary_id"].map(cell_id_dict)

        # CASE B: CosMx (The centroids hold the Cell ID string, but they are generated 1:1)
        elif len(boundary_list) == len(cell_id_dict):
            logger.info(
                "DEBUG: Mapping Cell IDs via 1:1 order alignment (Centroids to Boundaries)."
            )
            # Sort the internal IDs to guarantee perfect 1:1 alignment
            sorted_bounds = sorted(boundary_list)
            sorted_cents = sorted(list(cell_id_dict.keys()))
            b_to_c = dict(zip(sorted_bounds, sorted_cents))

            df_morph["Cell_ID"] = df_morph["boundary_id"].map(
                lambda b: cell_id_dict.get(b_to_c.get(b))
            )
        else:
            logger.warning(
                "DEBUG: Mismatch between number of Cell IDs and Boundaries. Cannot map safely!"
            )
            df_morph["Cell_ID"] = None

        # 5. Extract all target metrics
        target_labels = [
            "Area (µm²)",
            "Perimeter (µm)",
            "Convexity",
            "Circularity",
            "Principle axis angle (rad)",
        ]

        for label in target_labels:
            if label in domain.labels:
                vals, idxs = ms.query.get_labels(domain, label)
                metric_dict = dict(zip(list(idxs), list(vals)))
                df_morph[label] = df_morph["boundary_id"].map(metric_dict)

        # 5.5 Extract the Cluster labels
        # The variable 'chosen_cluster' is passed into this function (e.g., 'leiden_n10_r0.1')
        if chosen_cluster in domain.labels:
            cluster_vals, cluster_idxs = ms.query.get_labels(domain, chosen_cluster)
            cluster_dict = dict(zip(list(cluster_idxs), list(cluster_vals)))
            df_morph["Cluster"] = df_morph["boundary_id"].map(cluster_dict)
        else:
            logger.warning(
                f"Cluster label '{chosen_cluster}' not found in domain. Cells will be marked 'Unknown'."
            )
            df_morph["Cluster"] = "Unknown"

        # 6. Clean up and Save
        df_morph = df_morph.dropna(subset=["Cell_ID"])
        df_morph = df_morph.drop(columns=["boundary_id"])

        csv_path = out_dir / "morphometrics.csv"
        df_morph.to_csv(csv_path, index=False)
        logger.info(f"Successfully exported morphometrics to {csv_path}")

    except Exception as e:
        logger.warning(f"Failed to export Morphometrics CSV: {e}")


if __name__ == "__main__":
    # CosMx
    # muspan_domain_path = "C:/Users/bunga/python/Project2/CosMx/Spatial-Transcriptomics-CosMx-Xenium/analysis_Kitam/6_MuSpan/muspan_object_dnaku.muspan"
    # domain = ms.io.load_domain(muspan_domain_path)
    # chosen_cluster = "leiden_n50_r1.0"
    # cluster_list = ["0", "1", "2", "3", "6", "14"]

    # run_muspan_shapes(domain, chosen_cluster, cluster_list=None)

    # Xenium
    muspan_domain_path = "C:/Users/bunga/python/Project2/CosMx/Spatial-Transcriptomics-CosMx-Xenium/analysis_Tisku/6_MuSpan/muspan_object_pomni.muspan"
    domain = ms.io.load_domain(muspan_domain_path)
    chosen_cluster = "leiden_n10_r0.1"
    cluster_list = ["0", "1", "2"]

    run_muspan_shapes(domain, chosen_cluster, cluster_list)
