# Interactive Explorer

The Interactive Explorer is the primary dashboard for exploring your spatial dataset.

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/interactive" target="_blank">
    <img src="svk/assets/interactive.png"
         alt="Interactive tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

## Interface & Controls {#hide-me}

* Explore the physical spatial distribution (tissue map) and transcriptomic UMAP/PCA embeddings side-by-side. Navigation, zooming, and panning are automatically synchronized between both views.
* Use the Color By dropdown in the top control bar to color cells by Clusters, Annotations, or custom Metadata.
* Allows cells to be colored by a continuous heatmap of that gene's expression across the tissue and UMAP.
* Use the lasso tool in either the spatial plot or the embedding view to draw a boundary around a group of cells.
* Click the camera icons at the top right to download PNGs of the UMAP or Spatial tissue maps with custom background colors and legends.

## Label Composition Pie Chart {#hide-me}
Located in the bottom-right corner, the composition pie chart provides a breakdown of the currently active annotation/cluster.

* The chart dynamically updates based on the Slide and Sample selected in the top filter bar.
* Click on any specific slice in the pie chart. This acts as a global filter, isolating and highlighting only that specific cell type in the spatial and UMAP plots above, fading all other cells to the background.