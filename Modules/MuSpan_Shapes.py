"""Muspan module - Shapes analysis."""

import logging
import warnings
from pathlib import Path
from typing import List, Optional, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

try:
    import muspan as ms
except ModuleNotFoundError as err:
    logger.error("Could not load MuSpAn.")
    raise err


def ms_shapedesc(domain: Any) -> Any:
    """Calculates shape descriptors and generates a 2x2 plot."""
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


def ms_PrAxis(
    domain: Any, chosen_cluster: str, cluster_list: Optional[List[str]] = None
) -> Any:
    """Calculates shape orientation via principle axis and generates density/spatial plots."""
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

    # Calculate angles
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

    # Set up the figure
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

    # Visualise domain with principle axis
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

    # Visualise the chosen_cluster
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


def _export_morphometrics_csv(domain: Any, out_dir: Path, chosen_cluster: str) -> None:
    """Extracts shape metrics and maps them to original string Cell IDs."""
    logger.info("Extracting shape metrics into a CSV for web visualization...")

    try:
        cell_ids, cell_id_obj_indices = ms.query.get_labels(domain, "Cell ID")
        cell_id_dict = dict(zip(cell_id_obj_indices, cell_ids))

        # Because Area is only calculated on boundaries, this perfectly isolates them
        areas, boundary_indices = ms.query.get_labels(domain, "Area (µm²)")
        boundary_list = list(boundary_indices)

        df_morph = pd.DataFrame({"boundary_id": boundary_list})

        # Map logic (CosMx vs Xenium architecture diff)
        if boundary_list and boundary_list[0] in cell_id_dict:
            df_morph["Cell_ID"] = df_morph["boundary_id"].map(cell_id_dict)
        elif len(boundary_list) == len(cell_id_dict):
            sorted_bounds = sorted(boundary_list)
            sorted_cents = sorted(list(cell_id_dict.keys()))
            b_to_c = dict(zip(sorted_bounds, sorted_cents))
            df_morph["Cell_ID"] = df_morph["boundary_id"].map(
                lambda b: cell_id_dict.get(b_to_c.get(b))
            )
        else:
            logger.warning(
                "Mismatch between Cell IDs and Boundaries. Cannot map safely!"
            )
            df_morph["Cell_ID"] = None

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
                df_morph[label] = df_morph["boundary_id"].map(dict(zip(idxs, vals)))

        if chosen_cluster in domain.labels:
            cluster_vals, cluster_idxs = ms.query.get_labels(domain, chosen_cluster)
            df_morph["Cluster"] = df_morph["boundary_id"].map(
                dict(zip(cluster_idxs, cluster_vals))
            )
        else:
            df_morph["Cluster"] = "Unknown"

        df_morph = df_morph.dropna(subset=["Cell_ID"]).drop(columns=["boundary_id"])

        csv_path = out_dir / "morphometrics.csv"
        df_morph.to_csv(csv_path, index=False)
        logger.info(f"Successfully exported morphometrics to {csv_path}")

    except Exception as e:
        logger.error(f"Failed to export Morphometrics CSV: {e}")


def run_muspan_shapes(
    module_dir: Path,
    domain: Any,
    chosen_cluster: str,
    selected_celltypes: Optional[List[str]] = None,
) -> None:
    """Main orchestrator for shape analysis."""
    out_dir = module_dir / domain.name
    shape_dir = out_dir / "shape_analysis"
    shape_dir.mkdir(parents=True, exist_ok=True)

    # 1. Shape descriptors
    try:
        plt_fig = ms_shapedesc(domain)
        plt_fig.tight_layout()
        plt_fig.savefig(shape_dir / "ms_shapes.png")
        plt_fig.close()
        logger.info("Shape descriptor plot saved.")
    except Exception as e:
        logger.warning(f"Failed to generate shape descriptor plot: {e}")

    # 2. Principle axis
    try:
        plt_fig_axis = ms_PrAxis(domain, chosen_cluster, selected_celltypes)
        plt_fig_axis.tight_layout()
        plt_fig_axis.savefig(shape_dir / "ms_praxis.png")
        plt_fig_axis.close()
        logger.info("Shape orientation plot saved.")
    except Exception as e:
        logger.warning(f"Failed to generate shape orientation plot: {e}")

    # 3. Export CSV
    _export_morphometrics_csv(domain, shape_dir, chosen_cluster)
