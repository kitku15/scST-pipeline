# Cell-Cell Communication (CellPhoneDB)

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/ccc" target="_blank">
    <img src="svk/assets/ccc.png"
         alt="ccc tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

The **Cell-Cell Communication** tab visualizes the complex network of ligand-receptor interactions occurring across different cell populations identified by [CellPhoneDB](https://cellphonedb.readthedocs.io/en/latest/RESULTS-DOCUMENTATION.html#method-2-statistical-inference-of-interaction-specificity).

## Chord Diagram {#hide-me}

The center of this module is a chord diagram representing directed interactions:

*   The base of a connecting ribbon represents the **Sender** (the cell expressing the Ligand), pointing toward the **Receiver** (the cell expressing the Receptor).
*   The width of the ribbon correlates with the strength/significance of the interaction.
*   Click the 'Export' button to download a high-resolution PNG of the chord plot with a legend for the participating cell types or interactions.

## Filtering Options {#hide-me}

To make sense of dense networks, use the top control panel to filter the data:

1.  **Microenvironment:** Restrict the analysis to interactions occurring strictly within a specific spatial niche / region of interest.
2.  **Focal Cell Type:** View only interactions where a specific cell type is acting as either the sender or receiver.
3.  **Ligand-Receptor Pairs:** Use the searchable multiselect dropdown to view only specific molecular interactions.
4.  **Min Cells:** Filter out noise by requiring a minimum number of participating cells before an interaction is rendered.

## Visual Settings {#hide-me}

Use the Color By dropdown to change how the ribbons are painted:

*   **Receiver:** trace incoming signals to a specific cell.
*   **Sender:** trace outgoing signals from a specific cell.
*   **Interaction:** Color ribbons by specific Ligand-Receptor pairs.