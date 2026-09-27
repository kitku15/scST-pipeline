# Transcription Factor Analysis

The **Transcription Factor (TF) Enrichment** tab visualizes regulatory network activities inferred by [DecoupleR](https://decoupler-py.readthedocs.io/en/latest/).

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/tf" target="_blank">
    <img src="svk/assets/tf.png"
         alt="tf tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

The right-hand panel displays a clustered heatmap of Z-scaled Transcription Factor activity scores per cell type.

*   Click the cell type buttons above the heatmap to filter rows.
*   Red indicates high predicted TF activity; Blue indicates low activity.

The left-hand panel contains an interactive Vitessce viewer linked directly to the TF activity matrix.

*   Toggle between UMAP and Spatial views using the buttons at the top left.
*   Select a specific TF from the feature list (far right of the Vitessce panel) to color its predicted activity score onto the cells.
*   The leftmost viewer in this panel shows the standard categorical cell type labels for context, while the adjacent viewer shows the continuous TF activity.