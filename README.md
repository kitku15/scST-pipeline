# General Single-cell Spatial Transcriptomics Pipeline

An end-to-end Python pipeline for the processing, analysis, and visualization of Spatial Transcriptomics data. Built specifically to handle **NanoString CosMx** and **10x Genomics Xenium** datasets, this pipeline uses modern spatial data frameworks (spatialdata, scanpy, squidpy) and advanced spatial statistics (muspan).

## Key Features

-  Run the whole pipeline end-to-end, or execute specific modules independently.
- Supports both Machine Learning-based  ([CellTypist](https://www.celltypist.org/)) and Marker-based ([ScType](https://github.com/kris-nader/sc-type-py)) cell type annotation.
- Generates Delaunay, KNN, and Proximity graphs, along with cross-Pair Correlation Functions (PCF) via the [MuSpAn](https://www.muspan.co.uk/) library.
- Reproducible runs controlled entirely via a single TOML configuration file.
- Run on HPC with Singularity/Apptainer! (works nicely :D)
- Run locally using Docker! (if you must)
- Will integrate Nextflow soon (hopefully)
- Added runtime tracker to measure runtime each module

## Additional
This 22-week project was completed as part of the MRes in Bioinformatics and Theoretical Systems Biology at Imperial College London. It builds upon and adapts the [The ReCoDe-spatial-transcriptomics repository](https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics). The work was supervised by Tamas Korcsmaros and Balazs Bohar.

- Student: Bunga Tiasyaira Hutasuhut (Syaii) 
- Email: bth22@ic.ac.uk

