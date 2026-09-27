# Cell Type Annotation

The **Cell Type Annotation** tab allows you to track how cells are classified across different metadata columns, clustering resolutions, or annotations.

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/annotation" target="_blank">
    <img src="svk/assets/annotation.png"
         alt="annotation tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

At the center of this tab is an interactive Sankey diagram. It visualizes the flow of cells from one categorical state to another.

*   Use the dropdowns at the top to select the metadata columns you want to compare (e.g., Step 1: *Slide ID* → Step 2: *Cell Lineage* → Step 3: *Cell Type*).
*   Click the **+ Add Flow Step** button to introduce more layers to the comparison. You can remove steps using the **✕** button next to the respective dropdown.

The Insights panel on the right provides precise numerical breakdowns as you interact with the diagram:

*   Hovering over any cluster block reveals the total number of cells assigned to that label.
*   See where the cells in the current cluster originated from in the previous step, complete with cell counts and percentages.
*   See how a specific cluster splits apart and maps to new labels in the subsequent step.