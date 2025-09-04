

# SpaCNA: CNA detection from spatial transcriptomics using SpaCNA

**SpaCNA** is a computational pipeline for detecting Copy Number Alterations (CNAs) from spatial transcriptomics data by integrating three data modalities: **histology images**, **gene expression matrices**, and **spatial coordinates**.

This pipeline not only identifies CNAs but also includes downstream analysis modules for **estimating tumor content** and **detecting tumor boundaries**, providing a powerful toolkit for spatially-resolved studies of the tumor microenvironment.



-----

Of course. Here is the revised README section, updated in English with the specific version numbers you provided earlier.

-----

## Requirements & Installation

SpaCNA requires both **R** and **Python** environments. It is highly recommended to use a dedicated virtual environment (e.g., using `conda` or `venv`) to avoid conflicts with existing packages.

### Python Environment

  - **Python:** 3.7.12
  - **Key Packages:**
      - `numpy`==1.21.6
      - `matplotlib`==3.5.3
      - `torch`==1.12.1
      - `torchvision`==0.13.1
      - `pandas`==1.3.5
      - `scikit-learn`==1.0.2

It is recommended to install the required packages using the provided `requirements.txt` file.

**`requirements.txt`:**

```txt
numpy==1.21.6
matplotlib==3.5.3
torch==1.12.1
torchvision==0.13.1
pandas==1.3.5
scikit-learn==1.0.2
```

Install all dependencies with:

```bash
pip install -r requirements.txt
```

### R Environment

  - **R:** 4.2.1
  - **Key Packages:**
      - `Seurat` (4.2.0)
      - `biomaRt` (2.52.0)
      - `ComplexHeatmap` (2.12.1)
      - `parallelDist` (0.2.6)
      - `irlba` (2.3.5.1)
      - `edgeR` (3.38.4)
      - `ggplot2` (3.4.1)
      - `rootSolve` (1.8.2.3)
      - `patchwork` (1.1.2)
      - `glmnet` (4.1.8)

You can install them by running the following commands in your R console. This script handles packages from both CRAN and Bioconductor.

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install packages from CRAN
# For specific versions, you might need the 'remotes' package, e.g.:
# remotes::install_version("Seurat", version = "4.2.0")
install.packages(c("Seurat", "parallelDist", "irlba", "ggplot2", "rootSolve", "patchwork", "glmnet"))

# Install packages from Bioconductor
BiocManager::install(c("biomaRt", "ComplexHeatmap", "edgeR"))
```
-----
##  Usage Guide

The SpaCNA analysis workflow consists of the following main steps:

### Step 0: Data Preparation

Before you begin, please organize your files according to the following structure:

1.  **For Image Feature Extraction (Python):**
    Your sample directory should contain:

      * `tissue_hires_image.png`: The high-resolution histology image.
      * `spot.txt`: A single-column text file containing the cell barcodes for each spot.
      * `exp_location.txt`: A two-column text file (x, y) containing the pixel coordinates of each spot on the `tissue_hires_image.png`.



2.  **For CNA Analysis (R):**
    Your sample directory should contain:

    * `seurat_object.rds`: A standard Seurat object with the following structure:
        * The raw **gene expression matrix** must be stored at `obj@assays$Spatial@counts`.
        * It must contain **cell clustering results** (e.g., in `obj$seurat_clusters`) or custom cell-type annotations that can be used to define which spots will serve as the normal reference.

    * `gene_pos.rds`: A `data.frame` containing the genomic annotation for genes, which must include the following columns:
        * `symbol`: The gene symbol.
        * `chr`: The chromosome name (e.g., "chr1").
        * `start`: The gene's start position in base pairs (bp).
        * `end`: The gene's end position in base pairs (bp).

### Step 1: Extract Image Features

This step uses a pre-trained ResNet50 model to extract morphological features from the image tile centered on each spot.

  * **Script:** `get_spatial_feature.py`

  * **How to run:**

    1.  Open and edit the `get_spatial_feature.py` file.
    2.  Modify the `sample_dir` variable in the `main` block to point to the directory containing your image and coordinate files (must end with `/`).
    3.  Adjust the cropping radius `r` as needed based on image resolution.

    <!-- end list -->

    ```python
    if __name__ == "__main__":
        # Modify this to your sample directory
        sample_dir = "/your sample directory/"  
        
        # Cropping radius for each spot, adjustable depending on the resolution
        r = 50   

        # Run feature extraction
        get_pca_feature(sample_dir, r)
    ```

    4.  Execute the script from your terminal: `python get_spatial_feature.py`

  * **Output:**
    The script will automatically generate a `resnet50/` folder within your `sample_dir`, containing the extracted feature files, such as `features_pca.txt`.

### Step 2: Run SpaCNA for CNA Detection

This core step integrates gene expression, spatial coordinates, and image features to infer the CNA state for each spot.

  * **Script:** `SpaCNA.R`

  * **How to run:**
    In your R environment, load the `SpaCNA` function and run it with the following parameters.

    ```r
    # Example run
    sample_dir <- "/your sample directory/"
    image_dir <- "/your image feature directory/"
    plot_dir <- "/your to save results directory/"
    normal_clusters = c(cluster1, cluster2)

    cna_list <- SpaCNA(
        sample_dir = sample_dir,         # Directory containing seurat_object.rds and gene_pos.rds
        image_dir = image_dir,           # Directory containing the resnet50/ folder from Step 1
        plot_dir = plot_dir,             # Directory to save results and plots
        normal_clusters = normal_clusters # normal_clusters in `seurat_object.rds` that represent normal cells to be used as reference.
        # ... other advanced parameters
    )
    ```
    * **Output:**
    The script will generate the CNA results and plots within your `image_dir`.


### Step 3: Estimate Tumor Content (Optional)

This downstream analysis module estimates the proportion of tumor cells within each spot.

  * **Script:** `estimate_tumor_content.R`

  * **How to run:**

    ```r
    # Estimate tumor content
    seurat_obj_updated <- estimate_tumor_content(
        sample_dir = "/path/to/your/sample/data/",           # Directory containing seurat_object.rds
        spacna_dir = "/path/to/your/spacna/results/",        # Directory containing SpaCNA results (e.g., cns.rds)
        plot_dir = "/path/to/your/plot/output/",             # Directory for plot outputs
        K = 7                                               # Number of clones for clustering
    )
    ```

  * **Output:**
    Returns an updated Seurat object with a `tumor_content` column added to its `meta.data`.

### Step 4: Detect Tumor Edge (Optional)

This module identifies the boundary between tumor regions and normal tissue.

  * **Script:** `estimate_tumor_edge.R`

  * **How to run:**

    ```r
    # Detect tumor edge
    seurat_obj_final <- estimate_tumor_edge(
        sample_dir = "/path/to/your/sample/data/",
        spacna_dir = "/path/to/your/spacna/results/",
        plot_dir = "/path/to/your/plot/output/",
        tumor_content_dir = "/path/to/your/tumor_content/results/", # Directory with tumor content results
        # ... other parameters
    )
    ```

  * **Output:**
    Returns an updated Seurat object with `edge_score` and `edge` columns added to its `meta.data`.

## Demo Data

The `demo_data/` folder contains a sample dataset and the corresponding expected output files.

---

**Note:** The code and demo data were tested on a standard x64 Windows PC with an Intel CPU and 8 GB of RAM. Under these conditions, the runtime is ~10 minutes.