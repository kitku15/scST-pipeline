"""Muspan module - Cell-Cell Proximity analysis."""

import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
logger = getLogger(__name__)

try:
    import muspan as ms
except ModuleNotFoundError as err:
    logger.info("Could not load MuSpAn.")
    raise err


def create_proximity_network(
    domain, collection_name, network_name="proximity_network", max_edge_distance=1.0
):
    logger.info(
        f"Generating proximity network: '{network_name}' for collection: '{collection_name}'..."
    )
    qCells = ms.query.query(domain, ("Collection",), "is", collection_name)
    ms.networks.generate_network(
        domain,
        network_name=network_name,
        objects_as_nodes=qCells,
        network_type="Proximity",
        min_edge_distance=0,
        max_edge_distance=max_edge_distance,
    )


def plot_zoomed_khop_neighborhood(
    domain,
    network_name,
    collection_name,
    label_name,
    point_of_interest_id,
    max_distance,
    k=1,
    zoom_radius=30,
):
    khop_neighborhoods = ms.networks.khop_neighbourhood(
        domain, network_name=network_name, source_objects=[point_of_interest_id], k=k
    )
    neighborhood_nodes = list(khop_neighborhoods[point_of_interest_id])
    object_positions = {
        v: domain.objects[v].centroid for v in list(domain.objects.keys())
    }

    fig, ax = plt.subplots(figsize=(6, 6))
    qCells = ms.query.query(domain, ("Collection",), "is", collection_name)

    ms.visualise.visualise(
        domain,
        objects_to_plot=qCells,
        show_boundary=False,
        marker_size=50,
        add_cbar=False,
        ax=ax,
        shape_kwargs={"alpha": 0.2, "linewidth": 0.75, "facecolor": [0.7, 0.7, 0.7, 1]},
        scatter_kwargs={
            "linewidth": 0.75,
            "edgecolor": "black",
            "c": [0.7, 0.7, 0.7, 1],
        },
    )

    ms.visualise.visualise(
        domain,
        color_by=("label", label_name),
        objects_to_plot=neighborhood_nodes,
        show_boundary=False,
        marker_size=100,
        add_cbar=False,
        ax=ax,
        shape_kwargs={"alpha": 0.6, "linewidth": 0.75},
        scatter_kwargs={"linewidth": 0.75, "edgecolor": "black"},
    )

    ms.visualise.visualise(
        domain,
        color_by=("label", label_name),
        objects_to_plot=[point_of_interest_id],
        show_boundary=False,
        marker_size=150,
        add_cbar=False,
        ax=ax,
        shape_kwargs={"alpha": 1.0, "linewidth": 1.5},
        scatter_kwargs={"linewidth": 1.5, "edgecolor": "black"},
    )

    edge_list_me = [(point_of_interest_id, v) for v in neighborhood_nodes]
    nx.draw_networkx_edges(
        domain.networks[network_name],
        object_positions,
        edgelist=edge_list_me,
        edge_color="black",
        width=2,
        ax=ax,
        arrows=False,
    )

    this_point = domain.objects[point_of_interest_id].centroid
    ax.set_xlim(this_point[0] - zoom_radius, this_point[0] + zoom_radius)
    ax.set_ylim(this_point[1] - zoom_radius, this_point[1] + zoom_radius)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.title(
        f"Cell ID {point_of_interest_id} - {k}-hop Neighborhood - d={max_distance}"
    )
    plt.tight_layout()
    return plt


def calculate_contact_statistics(
    domain, network_name, label_name, clusters_of_interest
):
    logger.info("Calculating neighborhood contacts and compositions...")
    labels_array, object_indices = ms.query.get_labels(domain, label_name)
    label_categories = domain.labels[label_name]["categories"]

    degree_dict = {}
    avg_contacts = np.zeros((len(label_categories), len(clusters_of_interest) + 1))

    for idex, cat in enumerate(label_categories):
        this_id_query = ms.query.query(domain, ("label", label_name), "is", cat)
        khop_neighborhoods = ms.networks.khop_neighbourhood(
            domain, network_name=network_name, source_objects=this_id_query, k=1
        )
        degree_dict[cat] = [len(khop_neighborhoods[k]) for k in khop_neighborhoods]

        if len(khop_neighborhoods) > 0:
            these_compositions = np.zeros(
                (len(khop_neighborhoods), len(clusters_of_interest) + 1)
            )
            for kid, k in enumerate(khop_neighborhoods):
                neighbors = khop_neighborhoods[k]
                for n in neighbors:
                    idx_in_array = np.where(object_indices == n)[0][0]
                    neighbor_label = labels_array[idx_in_array]
                    if neighbor_label in clusters_of_interest:
                        col_idx = clusters_of_interest.index(neighbor_label)
                        these_compositions[kid, col_idx] += 1
                    else:
                        these_compositions[kid, -1] += 1

                total_contacts = np.sum(these_compositions[kid, :])
                if total_contacts > 0:
                    these_compositions[kid, :] = (
                        these_compositions[kid, :] / total_contacts
                    )
            avg_contacts[idex, :] = np.mean(these_compositions, axis=0)

    return degree_dict, avg_contacts, label_categories


def plot_contact_statistics(
    degree_dict,
    avg_contacts,
    label_categories,
    clusters_of_interest,
    max_distance,
    cluster_colors=None,
):
    if len(label_categories) <= 12:
        box_colors = sns.color_palette("Set3", len(label_categories))
    else:
        box_colors = sns.color_palette("husl", len(label_categories))

    if cluster_colors is None:
        if len(clusters_of_interest) <= 12:
            cluster_colors = sns.color_palette("Set3", len(clusters_of_interest))
        else:
            cluster_colors = sns.color_palette("husl", len(clusters_of_interest))

    bar_colors = list(cluster_colors) + [[0.8, 0.8, 0.8, 1]]
    fig, ax = plt.subplots(
        figsize=(18, 6), nrows=1, ncols=2, gridspec_kw={"width_ratios": [1, 2]}
    )
    data_to_plot = [degree_dict.get(v, []) for v in label_categories]

    sns.boxplot(
        data=data_to_plot,
        ax=ax[0],
        palette=box_colors,
        showfliers=False,
        saturation=0.75,
        boxprops={"alpha": 0.6},
        orient="v",
    )
    sns.stripplot(
        data=data_to_plot,
        ax=ax[0],
        palette=box_colors,
        jitter=0.3,
        alpha=0.3,
        orient="v",
    )

    ax[0].set_xticks(range(len(label_categories)))
    ax[0].set_xticklabels(label_categories, rotation=45, ha="right", fontsize=9)
    ax[0].set_ylabel("Number of Contacts (Degree)")
    ax[0].set_title("Total Contacts per Cell Type")

    for i, cluster in enumerate(clusters_of_interest):
        ax[1].bar(
            label_categories,
            avg_contacts[:, i],
            label=cluster,
            bottom=np.sum(avg_contacts[:, :i], axis=1),
            color=bar_colors[i],
        )
    ax[1].bar(
        label_categories,
        avg_contacts[:, -1],
        label="Other",
        bottom=np.sum(avg_contacts[:, :-1], axis=1),
        color=bar_colors[-1],
    )

    ax[1].set_xlabel("Cell Type")
    ax[1].set_ylabel("Proportion of Contacts")
    ax[1].set_title(f"Average Neighborhood Composition (d = {max_distance})")
    ax[1].set_xticks(range(len(label_categories)))
    ax[1].set_xticklabels(label_categories, rotation=45, ha="right", fontsize=9)
    ax[1].legend(title="Neighbor Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()
    return plt


def plot_global_network(
    domain,
    network_name,
    collection_name,
    centroid_collection_name,
    label_name,
    max_distance,
    figsize=(10, 8),
    marker_size=7.5,
):
    logger.info(f"Plotting global network '{network_name}'...")
    qCells = ms.query.query(domain, ("Collection",), "is", collection_name)
    qPoints = ms.query.query(domain, ("Collection",), "is", centroid_collection_name)
    points_and_cells = ms.query.query_container()
    points_and_cells.add_query(qPoints, "OR", qCells)

    fig, ax = plt.subplots(figsize=figsize, nrows=1, ncols=1)
    ms.visualise.visualise_network(
        domain,
        network_name=network_name,
        edge_weight_name=None,
        edge_cmap="Greys",
        edge_width=1,
        edge_vmin=0,
        edge_vmax=0.5,
        add_cbar=False,
        ax=ax,
        visualise_kwargs={
            "color_by": ("label", label_name),
            "objects_to_plot": points_and_cells,
            "add_cbar": False,
            "marker_size": marker_size,
            "shape_kwargs": {"alpha": 0.2, "linewidth": 0.75},
            "scatter_kwargs": {"linewidth": 0.2, "edgecolor": "black"},
        },
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    plt.title(f"Global Proximity Network ({label_name}), d={max_distance}")
    plt.tight_layout()
    return plt


def run_ms_cellproximity(
    module_dir,
    domain,
    chosen_cluster,
    selected_celltypes,
    selection_name,
    max_distance=200,
    cellboundary_label="Cell boundaries",
    centroid_collection_name="Cell centroids",
    color_dict=None,
):
    out_dir = Path(module_dir) / selection_name
    out_dir.mkdir(parents=True, exist_ok=True)

    if chosen_cluster in domain.labels:
        all_labels_in_domain = list(domain.labels[chosen_cluster]["categories"])
    else:
        raise ValueError(f"Label '{chosen_cluster}' not found in domain.labels.")

    if not selected_celltypes:
        selected_celltypes = all_labels_in_domain
        logger.info(
            f"No specific cell types provided. Showing composition broken down by all {len(selected_celltypes)} cell types."
        )
    else:
        logger.info(
            f"Showing composition specifically for {len(selected_celltypes)} selected cell types."
        )

    network_name = f"{selection_name}_proxnet{max_distance}"

    # 1. Generate Network
    create_proximity_network(
        domain=domain,
        collection_name=cellboundary_label,
        network_name=network_name,
        max_edge_distance=max_distance,
    )

    # 2. Compute statistics
    degrees, avg_contacts_matrix, all_labels = calculate_contact_statistics(
        domain=domain,
        network_name=network_name,
        label_name=chosen_cluster,
        clusters_of_interest=selected_celltypes,
    )

    # 3. Plot the boxplots and stacked bar charts
    my_colors = (
        [color_dict.get(ct, "gray") for ct in selected_celltypes]
        if color_dict
        else sns.color_palette("Set3", len(selected_celltypes))
    )
    plt_contact = plot_contact_statistics(
        degrees,
        avg_contacts_matrix,
        all_labels,
        selected_celltypes,
        max_distance,
        my_colors,
    )
    prox_dir = out_dir / "proximity_analysis"
    prox_dir.mkdir(parents=True, exist_ok=True)
    plt_contact.savefig(prox_dir / "contact_statistics.png")
    plt_contact.close()

    # 4. Plot a zoomed-in plot around a cell
    network_nodes = list(domain.networks[network_name].nodes)
    sample_cell_query = ms.query.query(
        domain, ("label", chosen_cluster), "is", selected_celltypes[0]
    )
    label_ids = ms.query.interpret_query(sample_cell_query)
    valid_ids = list(set(network_nodes).intersection(set(label_ids)))

    if len(valid_ids) == 0:
        logger.warning(
            f"No cells with label '{selected_celltypes[0]}' found in '{network_name}'. Skipping zoom plot."
        )
    else:
        sample_cell_id = valid_ids[0]
        plt_zoom = plot_zoomed_khop_neighborhood(
            domain,
            network_name,
            cellboundary_label,
            chosen_cluster,
            sample_cell_id,
            max_distance,
            1,
            50,
        )
        plt_zoom.savefig(prox_dir / "zoomed_khop_neighborhood.png")
        plt_zoom.close()

    # 5. Plot the full spatial network colored by cluster labels
    plt_global = plot_global_network(
        domain,
        network_name,
        cellboundary_label,
        centroid_collection_name,
        chosen_cluster,
        max_distance,
        (12, 10),
        5.0,
    )
    plt_global.savefig(prox_dir / "global_proximity_network.png")
    plt_global.close()

    # 6. Export Proximity Matrices to CSV
    try:
        df_contact = pd.DataFrame(
            avg_contacts_matrix,
            index=all_labels,
            columns=selected_celltypes + ["Other"],
        )
        df_contact.index.name = "Source_CellType"
        df_contact.to_csv(prox_dir / "contact_composition.csv")
        logger.info("Exported contact_composition.csv")
    except Exception as e:
        logger.warning(f"Failed to export contact composition CSV: {e}")
