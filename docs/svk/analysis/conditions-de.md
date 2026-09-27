# Conditions DE Analysis

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/conditions-de" target="_blank">
    <img src="svk/assets/conditions-de.png"
         alt="conditions-de tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

The **Conditions DE** tab visualizes PyDEseq2 results, allowing you to perform pairwise comparisons (e.g., Healthy vs. Disease) *within* a specific cell population.

## Pairwise Volcano Plot {#hide-me}

1.  Select a **Cell Type**.
2.  Select a **Pairwise Comparison** (e.g., `Disease_vs_Healthy`).

Located next to the volcano plot, this table provides an immediate summary list of the top upregulated genes for the Test condition and the Reference condition.

## Violin Plots {#hide-me}

The bottom half of the screen contains three panels for visualizing gene distributions.

*   The violins are split: one side represents the Test condition, the other represents the Reference condition.
*   Dotted Lines indicate the mean expression value for each respective condition.
*   Toggle the checkbox at the top right to drop cells with 0 expression. This allows you to compare transcriptomic shifts only within the subset of cells that are actively expressing the target gene.