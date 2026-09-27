# Cell Type DE Analysis

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/de-analysis" target="_blank">
    <img src="svk/assets/de-analysis.png"
         alt="de-analysis tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

The **Differential Expression (DE) Analysis** tab helps you identify the marker genes identified in specific cell populations / clusters.

## Volcano Plot {#hide-me}

*   Use the top dropdown to select an annotation grouping, then choose a target cluster.
*   The volcano plot visualizes statistical significance (-Log10 P-Value) against Fold Change. The **top-right quadrant** (red dots) contains genes that are significantly upregulated in your selected cluster compared to all other cells.

## Marker Table {#hide-me}

The bottom-right panel provides a quick-reference table displaying the top defining genes for *every* cluster in the currently selected annotation set. The currently selected cluster is highlighted.

## Gene Violin Plot {#hide-me}

The top-right panel features a split violin plot. 

*   Use the **Gene 1** searchable dropdown to pick any gene. 
*   The plot will display the expression distribution of that gene across all clusters side-by-side, allowing you to verify the specificity of a marker gene.