

# SpaCNA: A spatial-aware framework for detecting copy number alterations from spatial transcriptomics

**SpaCNA** is a computational pipeline for detecting Copy Number Alterations (CNAs) from spatial transcriptomics data by integrating three data modalities: **histology images**, **gene expression matrices**, and **spatial coordinates**.

This pipeline not only identifies CNAs but also includes downstream analysis modules for **estimating tumor content** and **detecting tumor boundaries**, providing a powerful toolkit for spatially-resolved studies of the tumor microenvironment.


## Requirements & Installation

SpaCNA requires both **R** and **Python** environments. It is highly recommended to use a dedicated virtual environment (e.g., using `conda` or `venv`) to avoid conflicts with existing packages.

### Python Environment

  - **Python:** 3.13.11
  - **Key Packages:**
      - `torch`==2.6.0+cu124
      - `torchvision`==0.21.0+cu124
      - `opencv-python`==4.12.0.88
      - `numpy`==2.2.6
      - `matplotlib`==3.10.8
      - `pandas`==2.3.3
      - `scikit-learn`==1.8.0
      - `scipy`==1.16.3

It is recommended to install the required packages using the provided `requirements.txt` file.

**`requirements.txt`:**

```txt
torch==2.6.0+cu124
torchvision==0.21.0+cu124
--extra-index-url https://download.pytorch.org/whl/cu124
opencv-python==4.12.0.88
numpy==2.2.6
pandas==2.3.3
scikit-learn==1.8.0
matplotlib==3.10.8
scipy==1.16.3
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
      - `dlm` (1.1.6)

You can install them by running the following commands in your R console. This script handles packages from both CRAN and Bioconductor.

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install packages from CRAN
# For specific versions, you might need the 'remotes' package, e.g.:
# remotes::install_version("Seurat", version = "4.2.0")
install.packages(c("Seurat", "parallelDist", "irlba", "ggplot2", "rootSolve", "patchwork", "glmnet","dlm"))

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

    Execute the script from your terminal and adjust the cropping radius `radius` as needed based on image resolution.

    <!-- end list -->

    ```terminal
    python get_spatial_feature.py --sample_dir "/your sample directory/" --radius 50
    ```

  * **Output:**
    The script will automatically generate a `resnet50/` folder within your `sample_dir`, containing the extracted feature files, such as `features_pca.txt`. Additionally, a `plot/` folder will be created within your `sample_dir`, which contains graph images corresponding to different image thresholds.

### Step 2: Run SpaCNA for CNA Detection

This core step integrates gene expression, spatial coordinates, and image features to infer the CNA state for each spot.

  * **Script:** `SpaCNA.R`

  * **How to run:**
    In your R environment, execute the script from your terminal.

    ```terminal

    Rscript SpaCNA.R --sample_dir="/your/sample/directory" --image_dir="/your/image feature/directory" --plot_dir="/path/to/output" --normal_clusters="1,2"
    ```
    * **Output:**
    The script will generate the CNA results and plots within your `image_dir`.


### Step 3: Estimate Tumor Content (Optional)

This downstream analysis module estimates the proportion of tumor cells within each spot.

  * **Script:** `estimate_tumor_content.R`

  * **How to run:**

    ```terminal
    Rscript estimate_tumor_content.R --sample_dir="/path/to/your/sample/data" --plot_dir="/path/to/output" --spacna_dir="/path/to/your/spacna/results"
    ```

  * **Output:**
    Returns an updated Seurat object with a `tumor_content` column added to its `meta.data`.

### Step 4: Detect Tumor Edge (Optional)

This module identifies the boundary between tumor regions and normal tissue.

  * **Script:** `estimate_tumor_edge.R`

  * **How to run:**

    ```terminal
    # Detect tumor edge
    Rscript estimate_tumor_edge.R --sample_dir="/path/to/data" --plot_dir="/path/to/output" --spacna_dir="/path/to/your/spacna/results" --tumor_content_dir="/path/to/your/tumor_content/results"

    ```

  * **Output:**
    Returns an updated Seurat object with `edge_score` and `edge` columns added to its `meta.data`.

## Demo Data

The `demo_data/` folder contains a sample dataset and the corresponding expected output files.

---

**Note:** The code and demo data were tested on a standard x64 Windows PC with an Intel CPU and 8 GB of RAM. Under these conditions, the runtime is ~10 minutes.