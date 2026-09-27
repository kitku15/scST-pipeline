# Composition Analysis

The **Composition Analysis** tab dynamically builds stacked bar charts to help you quantify and compare tissue heterogeneity across different samples, clusters, or disease states.

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/composition" target="_blank">
    <img src="svk/assets/composition.png"
         alt="composition tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

You can build complex cross-tabulations using three simple steps located in the top control panel:

1.  **Filter By (Optional - Col A):** Restrict the analysis to a specific subset of cells. For example, you can choose to only analyze the composition of cells where `Cell Lineage` is equal to `Immune`.
2.  **X-Axis Group (Col C):** Defines the grouping on the X-axis. For example, selecting `Disease Type` or `Sample ID` will create a separate bar for each unique condition/sample.
3.  **Breakdown / Colors (Col B):** Defines what the colored slices inside each bar represent. For example, selecting `Cell Type` will color the bar based on the proportion of different cell types.

All stacked bars are normalized to **100% (Fraction = 1.0)**. This allows for immediate visual comparison of proportions, regardless of variations in total cell counts between different samples or slides.