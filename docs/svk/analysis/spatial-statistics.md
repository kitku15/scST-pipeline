# Spatial Statistics

The **Spatial Statistics** tab provides insights into the physical organization of your tissue.

<div align="center">
  <a href="https://spatial-vis-kit.vercel.app/stats" target="_blank">
    <img src="svk/assets/stats.png"
         alt="stats tab screenshot"
         width="300">
  </a>
  <br>
  <small><em>Click on the image to try the live demo!</em></small>
</div>

1. **Neighborhoods**: Displays a heatmap showing which cell types are significantly colocalized or avoiding each other.
2. **Pair correlation function**: Pair Correlation Function (PCF) graphs analyze how the likelihood of finding two cell types together changes with distance. Select an Interaction Pair (e.g., Cluster A ↔ Cluster B). Spikes above the red dashed line (g(r) = 1) indicate the exact radius (in µm) where the two cell types strongly interact. # change wording for this?
3. **Morphology**: Compare physical cell properties across different clusters using violin plots. Use the dropdown to switch between metrics (e.g., *Area*, *Perimeter*). Click legend items to toggle specific clusters on/off.
4. **Autocorrelation**: Identifies highly structured spatial patterns through Moran's I. 
5. **Network Centrality**: Bar charts showing Degree and Closeness centrality for different clusters.

The right panel features a Vitessce tissue map viewer.