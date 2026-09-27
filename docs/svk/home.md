
<style>
  .md-content h1 {
    display: none;
  }
</style>

<p>
  <img src="../assets/svk_logo_hor.svg" alt="svk_logo_hor" width="500">
</p>

[**Spatial-VisKit**](https://github.com/korcsmarosgroup/Spatial-VisKit) is an interactive, web-based visualization tool designed for exploring single-cell spatial transcriptomics (ST) datasets.

[Get Started](getting-started.md){ .md-button .md-button--primary }
<a href="https://spatial-vis-kit.vercel.app/interactive" target="_blank" class="md-button">
  Live Demo
</a>
---

It was built specifically to explore advanced ST analysis results generated from
<a href="https://github.com/korcsmarosgroup/Spatial-Transcriptomics-CosMx-Xenium" target="_blank" rel="noopener noreferrer">scSpatial-Kit</a>,
while also providing a "Lite" mode allowing users to load and explore their own custom ST datasets easily.

*   **Tissue & UMAP Exploration:** Navigate physical tissue maps synchronized side-by-side with UMAP embeddings.
*   **Cell-Cell Communication:** Visualize ligand-receptor interactions globally across clusters or mapped directly onto physical tissue space.
*   **Spatial Statistics:** Analyze neighborhoods, distances, morphology, and spatial autocorrelation.
*   **Gene & TF Activity:** Map gene expression and DecoupleR transcription factor activity onto single cells.
*   **Differential Expression:** Explore marker genes and condition-specific expression shifts.

## Operating Modes {#hide-me}

Spatial-VisKit operates in two distinct modes depending on your data source:

| Mode | Description | Available Tabs |
|---|---|---|
| **Full Mode** | Designed for datasets processed by the scSpatial-Kit pipeline. | All 12 Tabs |
| **Lite Mode** | Designed for exploring your own dataset. | Interactive Explorer, Annotation, Composition, Multiplex |

