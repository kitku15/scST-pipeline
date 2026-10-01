# Spatial CCC (LIANA)

The **Spatial CCC** tab visualizes Cell-Cell Communication occurring at the single-cell level using [LIANA+](https://liana.readthedocs.io/en/stable/tutorials/notebooks/bivariate.html).

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/spatial-ccc" target="_blank">
    <img src="svk/assets/spatial-ccc.png"
         alt="spatial-ccc tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

Use the **Spatial Interaction Target** searchable dropdown to select a specific interaction footprint to map onto the tissue. Click the camera icons to download high-resolution PNGs of the interaction score, ligand expression, or receptor expression maps. A continuous color gradient legend is automatically appended.


*   **`LR_` Prefixes:** Specific Ligand-Receptor pairs.
*   **`CCC_` Prefixes:** NMF Communication Factors/Signatures.

When a specific Ligand-Receptor pair is selected, the Vitessce viewer dynamically configures a multi-panel layout:

1.  **Top Panel (Interaction Score):** Shows the combined spatial colocalization score.
2.  **Bottom Left (Ligand):** Shows the expression of the Sender's Ligand.
3.  **Bottom Right (Receptor):** Shows the expression of the Receiver's Receptor.

!!! tip "Intensity Threshold"
    Each of the three panels has an independent **Intensity Threshold** slider in the right-hand settings menu. You can independently filter out background noise for the ligand, receptor, and final score to isolate high-confidence interactions.