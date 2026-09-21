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
    logger.info(
        "Could not load MuSpAn. Install with:\n"
        "    pip install 'recode_st[muspan]' @ git+https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
    )
    raise err


def create_proximity_network(
    domain, collection_name, network_name="proximity_network", max_edge_distance=1.0
):
    """
    Generates a distance-weighted proximity network for the objects in your domain.

    Parameters:
    - domain: The MuSpAn domain object.
    - collection_name (str): The name of the collection containing your cells (e.g., 'Cell boundaries' or 'Cells').
    - network_name (str): The name to assign to the generated network.
    - max_edge_distance (float): The maximum distance between objects to be considered 'in contact' or connected.
    """
    print(
        f"Generating proximity network: '{network_name}' for collection: '{collection_name}'..."
    )

    # Query objects to build the network on
    qCells = ms.query.query(domain, ("Collection",), "is", collection_name)

    # Generate the network
    ms.networks.generate_network(
        domain,
        network_name=network_name,
        objects_as_nodes=qCells,
        network_type="Proximity",
        # distance_weighted=True,
        min_edge_distance=0,
        max_edge_distance=max_edge_distance,
    )
    print("Network generated successfully.")


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
    """
    Plots a zoomed-in view of a specific cell and its k-hop neighborhood, showing the network edges.

    Parameters:
    - domain: The MuSpAn domain object.
    - network_name (str): The name of your generated proximity network.
    - collection_name (str): The collection name of your cells.
    - label_name (str): The label used for cell types/clusters (e.g., 'Cluster ID').
    - point_of_interest_id (int): The unique MuSpAn ID of the central cell to visualize.
    - k (int): Number of hops for the neighborhood (1 = direct contacts).
    - zoom_radius (float): How far (in spatial units/microns) to show around the cell.
    """
    # Get the k-hop neighborhood for the specific point
    khop_neighborhoods = ms.networks.khop_neighbourhood(
        domain, network_name=network_name, source_objects=[point_of_interest_id], k=k
    )
    neighborhood_nodes = list(khop_neighborhoods[point_of_interest_id])

    # Extract positions for NetworkX drawing
    object_positions = {
        v: domain.objects[v].centroid for v in list(domain.objects.keys())
    }

    # plot
    fig, ax = plt.subplots(figsize=(6, 6))
    qCells = ms.query.query(domain, ("Collection",), "is", collection_name)

    # Background cells
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

    # Highlight neighborhood cells
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

    # Highlight the central cell of interest
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

    # Draw Network Edges
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

    # Zoom into the specific region
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
    """
    Calculates the degree (number of contacts) and the composition of those contacts.

    Parameters:
    - domain: The MuSpAn domain object.
    - network_name (str): The name of your proximity network.
    - label_name (str): The label used for cell types (e.g., 'cell_type').
    - clusters_of_interest (list): Specific cell types you want to focus on.

    Returns:
    - degree_dict: Dictionary containing contact numbers for each cluster.
    - avg_contacts: Matrix of contact proportions.
    - label_categories: List of all unique cell type labels in the domain.
    """
    print("Calculating neighborhood contacts and compositions...")

    # Get all categories and object indices
    labels_array, object_indices = ms.query.get_labels(domain, label_name)
    label_categories = domain.labels[label_name]["categories"]

    degree_dict = {}
    avg_contacts = np.zeros((len(label_categories), len(clusters_of_interest) + 1))

    for idex, cat in enumerate(label_categories):
        # Query cells of the current category
        this_id_query = ms.query.query(domain, ("label", label_name), "is", cat)

        # Get their 1-hop neighborhoods
        khop_neighborhoods = ms.networks.khop_neighbourhood(
            domain, network_name=network_name, source_objects=this_id_query, k=1
        )

        # Record degrees (number of contacts per cell)
        degree_dict[cat] = [len(khop_neighborhoods[k]) for k in khop_neighborhoods]

        # Calculate contact compositions
        if len(khop_neighborhoods) > 0:
            these_compositions = np.zeros(
                (len(khop_neighborhoods), len(clusters_of_interest) + 1)
            )

            for kid, k in enumerate(khop_neighborhoods):
                neighbors = khop_neighborhoods[k]
                for n in neighbors:
                    # Find what cluster this neighbor belongs to
                    idx_in_array = np.where(object_indices == n)[0][0]
                    neighbor_label = labels_array[idx_in_array]

                    if neighbor_label in clusters_of_interest:
                        col_idx = clusters_of_interest.index(neighbor_label)
                        these_compositions[kid, col_idx] += 1
                    else:
                        these_compositions[kid, -1] += 1  # 'Other' category

                # Normalize to proportions
                total_contacts = np.sum(these_compositions[kid, :])
                if total_contacts > 0:
                    these_compositions[kid, :] = (
                        these_compositions[kid, :] / total_contacts
                    )

            # Average composition for cells in this category
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
    """
    Generates a figure with a boxplot (contact degrees for ALL cells)
    and a stacked barplot (contact composition highlighting clusters_of_interest + 'Other').

    Parameters:
    - degree_dict (dict): From calculate_contact_statistics.
    - avg_contacts (ndarray): From calculate_contact_statistics.
    - label_categories (list): All cell types in the dataset.
    - clusters_of_interest (list): Cell types to highlight.
    - cluster_colors (list): Optional. Colors corresponding to clusters_of_interest.
    """
    # Palette for ALL cell types in the boxplot
    if len(label_categories) <= 12:
        box_colors = sns.color_palette("Set3", len(label_categories))
    else:
        box_colors = sns.color_palette("husl", len(label_categories))

    # Palette for the stacked bars (selected cell types + 'Other')
    if cluster_colors is None:
        if len(clusters_of_interest) <= 12:
            cluster_colors = sns.color_palette("Set3", len(clusters_of_interest))
        else:
            cluster_colors = sns.color_palette("husl", len(clusters_of_interest))

    # Add grey for the "Other" category at the end of the stack
    bar_colors = list(cluster_colors) + [[0.8, 0.8, 0.8, 1]]

    fig, ax = plt.subplots(
        figsize=(18, 6), nrows=1, ncols=2, gridspec_kw={"width_ratios": [1, 2]}
    )

    # Subplot 1: Boxplot of Degrees
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

    # Subplot 2: Stacked Bar Chart of Contact Composition
    for i, cluster in enumerate(clusters_of_interest):
        ax[1].bar(
            label_categories,
            avg_contacts[:, i],
            label=cluster,
            bottom=np.sum(avg_contacts[:, :i], axis=1),
            color=bar_colors[i],
        )

    # Plot 'Other'
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
    """
    Visualizes the full proximity network across the entire spatial domain.

    Parameters:
    - domain: The MuSpAn domain object.
    - network_name (str): The name of your generated proximity network.
    - collection_name (str): The name of the collection containing your cells
                             (must be the same one used to build the network).
    - centroid_collection_name (str): The name of the collection containing cell centroids.
    - label_name (str): The label used for coloring cell types/clusters.
    - figsize (tuple): Size of the output figure (default is 10x8).
    - marker_size (float): Size of the cell markers/points.
    """
    print(f"Plotting global network '{network_name}'...")

    # 1. Query the boundaries you used to build the network
    qCells = ms.query.query(domain, ("Collection",), "is", collection_name)

    # 2. Query the cell centers (centroids)
    qPoints = ms.query.query(domain, ("Collection",), "is", centroid_collection_name)

    # 3. Create a query container to combine both cell centers and boundaries
    points_and_cells = ms.query.query_container()
    points_and_cells.add_query(qPoints, "OR", qCells)

    fig, ax = plt.subplots(figsize=figsize, nrows=1, ncols=1)

    # Visualize the network
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
    # Clean up axes for spatial plotting
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

    # Ensure the chosen label exists in the domain
    if chosen_cluster in domain.labels:
        all_labels_in_domain = list(domain.labels[chosen_cluster]["categories"])
    else:
        raise ValueError(f"Label '{chosen_cluster}' not found in domain.labels.")

    # If no cell types are provided, we show composition broken down by ALL cell types
    if not selected_celltypes:
        selected_celltypes = all_labels_in_domain
        logger.info(
            f"No specific cell types provided. Showing composition broken down by all {len(selected_celltypes)} cell types."
        )
    else:
        logger.info(
            f"Showing composition specifically for {len(selected_celltypes)} selected cell types. The rest will be grouped as 'Other'."
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
    if color_dict is not None:
        my_colors = [color_dict.get(ct, "gray") for ct in selected_celltypes]
    else:
        if len(selected_celltypes) <= 12:
            my_colors = sns.color_palette("Set3", len(selected_celltypes))
        else:
            my_colors = sns.color_palette("husl", len(selected_celltypes))

    plt_contact = plot_contact_statistics(
        degree_dict=degrees,
        avg_contacts=avg_contacts_matrix,
        label_categories=all_labels,
        clusters_of_interest=selected_celltypes,
        max_distance=max_distance,
        cluster_colors=my_colors,
    )
    plt_contact.savefig(f"{out_dir}/contact_statistics.png")
    plt_contact.close()

    # 4. Plot a zoomed-in plot around a cell
    network_nodes = list(
        domain.networks[network_name].nodes
    )  # Get the list of all cell IDs
    sample_cell_query = ms.query.query(
        domain, ("label", chosen_cluster), "is", selected_celltypes[0]
    )  # Query all cell IDs that belong to your target cluster
    label_ids = ms.query.interpret_query(sample_cell_query)
    valid_ids = list(
        set(network_nodes).intersection(set(label_ids))
    )  # Find IDs that have the label AND are in the network

    if len(valid_ids) == 0:
        logger.warning(f"No cells with label '{selected_celltypes[0]}' found in '{network_name}'. Skipping zoom plot.")
        return # Exits the function without crashing

    sample_cell_id = valid_ids[0]  # Grab the first valid cell ID
    print(
        f"Successfully found valid cell ID: {sample_cell_id} for cluster {selected_celltypes[0]}"
    )

    plt_zoom = plot_zoomed_khop_neighborhood(
        domain=domain,
        network_name=network_name,
        collection_name=cellboundary_label,
        label_name=chosen_cluster,
        point_of_interest_id=sample_cell_id,
        max_distance=max_distance,
        k=1,
        zoom_radius=50,
    )
    plt_zoom.savefig(f"{out_dir}/zoomed_khop_neighborhood.png")
    plt_zoom.close()

    # 5. Plot the full spatial network colored by cluster labels
    plt_global = plot_global_network(
        domain=domain,
        network_name=network_name,
        collection_name=cellboundary_label,
        centroid_collection_name=centroid_collection_name,
        label_name=chosen_cluster,
        max_distance=max_distance,
        figsize=(12, 10),
        marker_size=5.0,
    )
    plt_global.savefig(f"{out_dir}/global_proximity_network.png")
    plt_global.close()

    # 6. Export Proximity Matrices to CSV for aggregation
    try:
        df_contact = pd.DataFrame(
            avg_contacts_matrix,
            index=all_labels,
            columns=selected_celltypes + ["Other"],
        )
        df_contact.index.name = "Source_CellType"
        df_contact.to_csv(out_dir / "contact_composition.csv")
        logger.info("Exported contact_composition.csv")
    except Exception as e:
        logger.warning(f"Failed to export contact composition CSV: {e}")
