# Getting Started

## Prerequisites

To run this application, you need Docker installed on your computer.

1. Download and install **[Docker Desktop](https://www.docker.com/products/docker-desktop/)**.
2. Make sure Docker Desktop is open and running in the background before proceeding.

---

## Step 1: Prepare Your Data

Place your dataset folder into the `public/` directory of the Spatial-VisKit application.

1. Locate your dataset folder. This might be an output from the pipeline or your own dataset in spatial data `.zarr` format.
2. Copy that entire folder into the `public/` folder of Spatial-VisKit.
3. Your folder structure should look like this:

```text
Spatial-VisKit/
├── public/
│   ├── my_data/
│   │   ├── dataset_config.json
│   │   ├── my_data.zarr
│   │   ├── my_data_tf.zarr        # (Pipeline output only)
│   │   └── aux_data/              # (Pipeline output only)
├── .env
├── docker-compose.yml
└── ...
```

---

## Step 2: Configure the App

There are two configuration files you need to set up:

### A. The `.env` File {#hide-me}
Open the file named `.env` (located in the main Spatial-VisKit folder) using any basic text editor. Update the values to match the dataset you want to view.

```env title=".env"
# 1. The exact name of your dataset folder inside the /public directory
ACTIVE_DATASET_FOLDER=my_data

# 2. What mode is the app in? ("full" or "lite")
VITE_APP_MODE=full

# 3. A title to be displayed in the app header
VITE_PROJECT_TITLE="CosMx SMI: My Data XYZ"

# 4. Do not change this
VITE_API_BASE_URL=http://localhost:8000
```

### B. The `dataset_config.json` File {#hide-me}
Inside your specific dataset folder (e.g., `public/my_data/`), there must be a file named `dataset_config.json`. This tells the visualizer how to read your AnnData/Zarr arrays. 

Select the tab below that matches your operating mode to see a configuration example:

=== "Full Mode (Pipeline Output)"
    ```json title="public/my_data/dataset_config.json"
    {
      "zarr_filename": "my_data_web.zarr", 
      "tf_zarr_filename": "my_data_tf_web.zarr", 
      "spatial_key": "global", 
      "slide_col": "slide_ID", 
      "sample_col": "sample_id", 
      "primary_annotation": "Final_Annotation", 
      "dynamic_annotations": [
        {"name": "Cell Clusters (Leiden)", "prefix": "leiden" }
      ],
      "extra_obs_sets": [
        { "name": "Disease Type", "path": "obs/DiseaseType" },
        { "name": "Treatment Response", "path": "obs/TreatmentResponse" },
        { "name": "FOV", "path": "obs/fov" },
        { "name": "Cell Type", "path": "obs/Final_Annotation" }
      ]
    }
    ```

=== "Lite Mode (Bring Your Own Data)"
    ```json title="public/my_data/dataset_config.json"
    {
      "zarr_filename": "my_data.zarr",
      "spatial_key": "global",
      "slide_col": "slide_id",
      "sample_col": "sample_id",
      "primary_annotation": "leiden_n30_r1.0",
      
      "available_embeddings": [
        { "name": "UMAP", "path": "obsm/X_umap" },
        { "name": "PCA", "path": "obsm/X_pca" }
      ],
      
      "dynamic_annotations": [],
      
      "extra_obs_sets": [
        { "name": "Disease Type", "path": "obs/DiseaseType" },
        { "name": "Cell Type", "path": "obs/CellType" }
      ]
    }
    ```

---

## Step 3: Run the Application

Once your data is in the `public/` folder and your configuration files are saved:

1. Open your computer's Terminal (Mac/Linux) or Command Prompt/PowerShell (Windows).
2. Navigate to the main `Spatial-VisKit` folder.
3. Run the following command to build and start the containers:
   ```bash
   docker-compose up --build
   ```
4. Wait for the setup to finish. It will download the necessary dependencies and launch the backend and frontend.
5. Open your web browser and navigate to: [http://localhost:5173](http://localhost:5173)

### Stopping the Application {#hide-me}

To stop the application, go back to your terminal window where Docker is running and press `Ctrl + C`.

!!! warning "Changing Datasets"
    If you ever change the `ACTIVE_DATASET_FOLDER` in your `.env` file to look at a different dataset, you **must** stop the application (`Ctrl + C`) and restart it by running `docker-compose up --build` again so Docker can mount the new folder.
