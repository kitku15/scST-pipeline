"""Muspan module - spatial statistics and graph analysis."""

import json
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def calculate_and_plot_cross_pcf(
    ms,
    domain,
    cluster_labels,
    module_dir,
    unique_clusters,
    cell_types=("14", "18"),
    max_R=200,
    annulus_step=5,
    annulus_width=25,
    visualise_output=True,
):
    """Calculates and plots the cross-PCF for two selected cell types."""
    cell_type_1, cell_type_2 = str(cell_types[0]), str(cell_types[1])

    if cell_type_1 not in unique_clusters or cell_type_2 not in unique_clusters:
        logger.warning(
            f"Cell types '{cell_type_1}' or '{cell_type_2}' not found in unique clusters. Skipping specific cross-PCF plot."
        )
        return

    logger.info(f"Calculating cross-PCF for {cell_type_1} and {cell_type_2}...")

    # Query populations
    centroids_query = ms.query.query(domain, ("collection",), "is", "Cell centroids")
    pop_A = ms.query.query_container(
        centroids_query,
        "AND",
        ms.query.query(domain, ("label", cluster_labels), "is", str(cell_type_1)),
    )
    pop_B = ms.query.query_container(
        centroids_query,
        "AND",
        ms.query.query(domain, ("label", cluster_labels), "is", str(cell_type_2)),
    )

    try:
        # Calculate cross-PCF
        r, PCF = ms.spatial_statistics.cross_pair_correlation_function(
            domain=domain,
            population_A=pop_A,
            population_B=pop_B,
            max_R=max_R,
            annulus_step=annulus_step,
            annulus_width=annulus_width,
            visualise_output=visualise_output,
        )

        domain_name = domain.name
        stat_dir = module_dir / domain_name / "spatial_stats"
        stat_dir.mkdir(parents=True, exist_ok=True)
        # Save PCF plot
        pcf_plot_path = stat_dir / f"cross_pair_correlation_function_{cell_type_1}_{cell_type_2}.png"
        plt.savefig(pcf_plot_path)
        logger.info(f"Cross-PCF plot saved at {pcf_plot_path}")

        # Visualize and save cell type points
        query_1_2 = ms.query.query(
            domain,
            ("label", cluster_labels),
            "in",
            [str(cell_type_1), str(cell_type_2)],
        )

        fig, ax = ms.visualise.visualise(
            domain, color_by=("label", cluster_labels), objects_to_plot=query_1_2
        )

        ax.set_title(f"{cell_type_1} vs {cell_type_2}", fontsize=15)
        ax.tick_params(axis="both", which="major", labelsize=10)
        ax.set_xlabel(f"{cell_type_1}", fontsize=15)
        ax.set_ylabel(f"{cell_type_2}", fontsize=15)

        vis_plot_path = stat_dir / f"visualize_{cell_type_1}_{cell_type_2}.png"
        plt.savefig(vis_plot_path)
        logger.info(f"Visualization saved at {vis_plot_path}")

    except Exception as e:
        logger.warning(
            f"Failed to calculate or plot cross-PCF for {cell_type_1} vs {cell_type_2}: {e}"
        )
        plt.close("all")


def calculate_pairwise_cross_pcf(
    ms,
    domain,
    module_dir,
    cluster_labels: str,
    unique_clusters: list[str],
):
    """Calculates and plots the cross-PCF for all unique unordered cell type pairs."""

    num_clusters = len(unique_clusters)

    if num_clusters < 2:
        logger.warning(
            "Less than 2 unique clusters found. Skipping pairwise cross-PCF matrix."
        )
        return

    # Initialize JSON dictionary
    cross_pcf_json = {"r_distances": [], "pairs": {}}

    # Only draw the massive matrix plot if we have a reasonable number of clusters (< 15)
    plot_matrix = num_clusters <= 15
    if plot_matrix:
        fig, axes = plt.subplots(num_clusters, num_clusters, figsize=(40, 40))
    else:
        logger.warning(
            f"Too many clusters ({num_clusters}) for matrix plot. Skipping plot, but saving JSON data."
        )

    for i, cluster_i in enumerate(unique_clusters):
        for j, cluster_j in enumerate(unique_clusters):
            if j <= i:
                if plot_matrix:
                    axes[i, j].axis("off")  # Optional: remove lower triangle plots
                continue

            centroids_query = ms.query.query(
                domain, ("collection",), "is", "Cell centroids"
            )

            pop_A = ms.query.query_container(
                centroids_query,
                "AND",
                ms.query.query(domain, ("label", cluster_labels), "is", cluster_i),
            )
            pop_B = ms.query.query_container(
                centroids_query,
                "AND",
                ms.query.query(domain, ("label", cluster_labels), "is", cluster_j),
            )
            logger.info(f"Calculating cross-PCF: {cluster_i} vs {cluster_j}")

            try:
                r, pcf = ms.spatial_statistics.cross_pair_correlation_function(
                    domain,
                    pop_A,
                    pop_B,
                    max_R=200,
                    annulus_step=5,
                    annulus_width=10,
                )

                # Only save 'r' once since it's the same for all pairs
                if len(cross_pcf_json["r_distances"]) == 0:
                    cross_pcf_json["r_distances"] = r.tolist()

                cross_pcf_json["pairs"][f"{cluster_i}|{cluster_j}"] = np.nan_to_num(
                    pcf, nan=0.0
                ).tolist()

                if plot_matrix:
                    ax = axes[i, j]
                    ax.plot(r, pcf)
                    ax.axhline(1, color="k", linestyle=":")
                    ax.tick_params(axis="both", which="major", labelsize=15)
                    ax.set_ylabel(f"$g_{{{cluster_i},{cluster_j}}}(r)$", fontsize=20)
                    ax.set_xlabel("$r$", fontsize=20)

            except Exception as e:
                logger.warning(
                    f"Failed to calculate cross-PCF for {cluster_i} vs {cluster_j}: {e}"
                )
                if plot_matrix:
                    axes[i, j].axis("off")

    domain_name = domain.name
    stat_dir = Path(module_dir) / domain_name / "spatial_stats"
    stat_dir.mkdir(parents=True, exist_ok=True)
    output_path = stat_dir / "cross_pair_correlation_function_all.png"

    if plot_matrix:
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close(fig)
        logger.info(f"Cross-PCF matrix plot saved at {output_path}")

    # Write out the JSON file (This always runs, regardless of plotting)
    json_output_path = stat_dir / "cross_pcf_all.json"
    with open(json_output_path, "w") as f:
        json.dump(cross_pcf_json, f)

    logger.info(f"Cross-PCF matrix plot saved at {output_path}")


def run_muspan_stats(module_dir, domain, cluster_labels, cell_types):
    """Run Muspan spatial statistics analysis on Xenium data."""
    try:
        import muspan as ms
    except ModuleNotFoundError as err:
        logger.info("Could not load MuSpAn.")
        raise err

    module_dir.mkdir(exist_ok=True)

    # Get cluster labels
    labels_array, _ = ms.query.get_labels(domain, cluster_labels)
    unique_clusters = np.unique(labels_array).astype(str).tolist()
    logger.info(f"Found {len(unique_clusters)} unique cell types.")

    # Run pairwise analysis for specific pair
    logger.info("Calculating pairwise cross-PCF for selected cell types...")
    calculate_and_plot_cross_pcf(
        ms=ms,
        domain=domain,
        cluster_labels=cluster_labels,
        module_dir=module_dir,
        unique_clusters=unique_clusters,
        cell_types=cell_types,
    )

    # Full pairwise matrix
    logger.info("Calculating pairwise cross-PCF for all cell types...")
    calculate_pairwise_cross_pcf(
        ms=ms,
        domain=domain,
        module_dir=module_dir,
        cluster_labels=cluster_labels,
        unique_clusters=unique_clusters,
    )


def aggregate_muspan_across_conditions(
    module_6_dir,
    sample_condition_dict,
    reference_condition=None,
    base_selection_name=None,
):
    """Aggregates Morphometrics and Proximity Networks across conditions."""
    import logging
    from pathlib import Path

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    logger = logging.getLogger(__name__)
    logger.info("=== Running MuSpAn Condition-Level Aggregations ===")

    agg_dir = Path(module_6_dir) / "Aggregated_Results"
    agg_dir.mkdir(exist_ok=True)

    all_morph = []
    contact_dict = {}

    # 1. Gather Data from all Sample Folders
    for sample, condition in sample_condition_dict.items():
        # Build the exact nested path using the selection name
        if base_selection_name:
            selection_name = f"{base_selection_name}_{sample}"
            sample_dir = Path(module_6_dir) / sample / selection_name
        else:
            sample_dir = Path(module_6_dir) / sample

        # Gather Morphometrics
        morph_file = sample_dir / "morphometrics.csv"
        if morph_file.exists():
            df = pd.read_csv(morph_file)
            df["Sample"] = sample
            df["Condition"] = condition
            all_morph.append(df)
        else:
            logger.warning(f"Could not find {morph_file}")

        # Gather Contacts
        contact_file = sample_dir / "proximity_analysis" / "contact_composition.csv"
        if contact_file.exists():
            df_cont = pd.read_csv(contact_file, index_col=0)
            if condition not in contact_dict:
                contact_dict[condition] = []
            contact_dict[condition].append(df_cont)
        else:
            logger.warning(f"Could not find {contact_file}")

    # 2. Plot Morphometrics (Area & Circularity Boxplots)
    if all_morph:
        logger.info("Plotting Morphometrics...")
        df_morph = pd.concat(all_morph)
        for metric in ["Area (µm²)", "Circularity"]:
            if metric in df_morph.columns:
                plt.figure(figsize=(14, 6))
                sns.boxplot(
                    data=df_morph,
                    x="Cluster",
                    y=metric,
                    hue="Condition",
                    palette="Set2",
                    showfliers=False,
                )
                plt.xticks(rotation=45, ha="right")
                plt.title(f"Cell {metric} across Conditions")
                plt.tight_layout()
                plt.savefig(agg_dir / f"Aggregated_{metric.split(' ')[0]}.png", dpi=300)
                plt.close()
        df_morph.to_csv(agg_dir / "All_Samples_Morphometrics.csv", index=False)

    # 3. Plot Consensus Proximity (Contact Composition Heatmaps)
    if contact_dict:
        logger.info("Plotting Contact Networks...")
        mean_contacts = {}
        for cond, df_list in contact_dict.items():
            # Average the matrices for this condition
            concat_df = pd.concat(df_list)
            mean_df = concat_df.groupby(concat_df.index).mean()
            mean_contacts[cond] = mean_df

            plt.figure(figsize=(10, 8))
            sns.heatmap(
                mean_df,
                annot=True,
                fmt=".2f",
                cmap="YlGnBu",
                cbar_kws={"label": "Proportion of Contacts"},
            )
            plt.title(f"Average Cell-Cell Contact Composition: {cond}")
            plt.ylabel("Source Cell")
            plt.xlabel("Neighbor Cell")
            plt.tight_layout()
            plt.savefig(agg_dir / f"Consensus_Contacts_{cond}.png", dpi=300)
            plt.close()

        # Plot Difference Heatmap using the user-defined reference condition!
        if reference_condition and reference_condition in mean_contacts:
            for cond, mean_df in mean_contacts.items():
                if cond != reference_condition:
                    diff_df = mean_df - mean_contacts[reference_condition]
                    plt.figure(figsize=(10, 8))
                    vmax = diff_df.abs().max().max()
                    sns.heatmap(
                        diff_df,
                        annot=True,
                        fmt=".2f",
                        cmap="coolwarm",
                        center=0,
                        vmin=-vmax,
                        vmax=vmax,
                    )
                    plt.title(f"Δ Contact Composition ({cond} - {reference_condition})")
                    plt.tight_layout()
                    plt.savefig(
                        agg_dir / f"Diff_Contacts_{cond}_vs_{reference_condition}.png",
                        dpi=300,
                    )
                    plt.close()

    logger.info(f"MuSpAn aggregations saved to {agg_dir}")
