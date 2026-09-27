# Multiplex Gene Overlay

The **Multiplex Overlay** tab transforms your spatial transcriptomics data into an interface resembling immunofluorescence (IF) microscopy. It allows you to visualize the spatial overlap of up to 5 specific genes simultaneously.

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/multiplex" target="_blank">
    <img src="svk/assets/multiplex.png"
         alt="multiplex tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>


!!! warning "Sample Selection Required"
    Because this view plots physical X/Y coordinates, you **must** select a specific *Sample* from the top dropdown. This prevents cells from different tissues from overlapping on top of one another.

*   Click **Add Channel** to introduce a new gene overlay. You can assign each gene a distinct primary color (Red, Green, Blue, Magenta, Cyan, Yellow).
*   When multiple genes are expressed in the same cell, their colors blend additively. For example, if a cell expresses both a Red-assigned gene and a Green-assigned gene, the resulting cell will appear Yellow.
*   Each channel features an *Intensity Threshold* slider. Adjust this to filter out low-expression background noise, ensuring that only cells with strong signals contribute to the color mapping.
*   Adjust the global point size slider to make cells larger and more visible, or smaller to distinguish densely packed regions.
*   You can toggle the visibility of "Background Cells" (cells not expressing any of your selected genes) to clean up the plot.