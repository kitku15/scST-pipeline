# Installation Guide

!!! abstract "Overview"
    scSpatial-Kit is fully containerized using **Docker**. You do not need to install Python, Conda, or any libraries directly on your host machine.

## 1. Prerequisites
You only need to have [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) installed on your machine or HPC cluster. 

*(When using singularity, you need to convert the Docker image to an [Apptainer/Singularity](https://apptainer.org/docs/user/main/docker_and_oci.html) image).*

## 2. Clone the Repository
Download the pipeline code to your local machine or server.

```bash
git clone https://github.com/korcsmarosgroup/Spatial-Transcriptomics-CosMx-Xenium.git
cd Spatial-Transcriptomics-CosMx-Xenium
```

## 3. Build the Docker Image
Build the Docker image using the provided `Dockerfile`. This process will download the base OS and install all required bioinformatics packages (Scanpy, Squidpy, CellPhoneDB, LIANA, MuSpAn, etc.). 

*(Note: This step may take 10-20 minutes depending on your internet connection, but you only have to do it once).*

```bash
docker build -t scSpatial-Kit:latest .
```

!!! info "GPU Acceleration (Optional)"
    If your machine has an NVIDIA GPU and you want to accelerate the deep learning steps (Module 2's scVI/scVIVA models), ensure you have the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) installed on your host machine.

## 4. Running the Container
Because the pipeline processes large amounts of data, you must mount your local directory into the container so it can read data and save outputs. The easiest way to use the pipeline is to launch an interactive shell inside the container:

```bash
docker run -it --rm \
  -v $(pwd):/pipeline \
  -w /pipeline \
  scSpatial-Kit:latest /bin/bash
```

* **`-it`**: Runs the container interactively.
* **`--rm`**: Deletes the container instance when you exit (keeps your system clean).
* **`-v $(pwd):/pipeline`**: Maps your current folder into the container at `/pipeline`. Any files the pipeline creates will appear on your actual computer.

## 5. Verify Installation
Once inside the Docker container's terminal, you can verify everything works by triggering a "dry run" using the template config file:

```bash
# You are now running commands inside the Docker container
python run_pipeline.py config_CosMx.toml --dry-run
```
If the image was built correctly, Snakemake will print the execution DAG without throwing any import errors. You are now ready to format your data and run the pipeline!
