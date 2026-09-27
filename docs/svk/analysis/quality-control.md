# Quality Control

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/qc" target="_blank">
    <img src="svk/assets/qc.png"
         alt="qc tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

The **QC tab** allows you to inspect the pre-filter quality metrics of your dataset. The dashboard displays histograms of raw cell metrics before filtering was applied:

1.  Total Transcripts per Cell
2.  Unique Genes per Cell 
3.  Cell Area
4.  Mean DAPI / Nucleus Ratio

These vertical lines indicate the exact cutoff thresholds that were applied during the pipeline's QC filtering step. Use the **Slide** dropdown at the top to switch between an aggregate view of all slides, or inspect individual slides to see if a specific batch had significantly lower quality metrics.