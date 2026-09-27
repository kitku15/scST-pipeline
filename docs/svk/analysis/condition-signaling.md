# Condition Signaling (Causal)

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/conditions-causal" target="_blank">
    <img src="svk/assets/conditions-causal.png"
         alt="conditions-causal tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

The **Causal Network** tab combines [LIANA+](https://liana.readthedocs.io/en/stable/tutorials/notebooks/targeted.html) and [Corneto](https://corneto.org/stable/index.html) outputs to map how external signals (Ligands/Receptors) trigger intracellular cascades to regulate Transcription Factors under specific disease or treatment conditions.

Use the top control bar to define the context of the signaling network:
1.  Select the condition contrast (e.g., treated vs. untreated).
2.  Select the cell population providing the Ligands.
3.  Select the cell population receiving the signal and undergoing intracellular changes.

### Condition-Altered Signals (Top-Left) {#hide-me}

A dot plot showing Ligand-Receptor pairs that are significantly shifted between your conditions.

*   **Red:** Up-regulated in the test condition.
*   **Blue:** Down-regulated in the test condition.
*   **Size:** Represents statistical significance.

### Altered TF Activity (Bottom-Left) {#hide-me}

A horizontal bar chart showing which Transcription Factors have statistically significant shifts in activity within the *Receiver Cell*.

### Causal Network (Right) {#hide-me}

An interactive, D3 force-directed graph mapping the prior knowledge pathways linking the active cell-surface Receptors to the affected nuclear Transcription Factors.

*   **Green Nodes:** Receptors
*   **Grey Nodes:** Intermediate Kinases / Intracellular proteins
*   **Purple Nodes:** Transcription Factors
*   **Solid Lines:** Biological Activation
*   **Dashed Lines:** Biological Inhibition

*You can click and drag nodes to explore the network structure.*