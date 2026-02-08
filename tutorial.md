# SpaCNA Demo: Step-by-Step Execution Guide

## System Requirements
- **Operating System:** Windows x64, macOS, or Linux
- **Memory:** Minimum 8 GB RAM
- **Storage:** 2 GB free space
- **Python:** 3.13.11
- **R:** 4.2.1

## 1. Environment Configuration

### 1.1 Python Environment Setup

Open your terminal/command prompt and execute:

```bash
# Create virtual environment
python -m venv spacna_demo_env

# Activate environment
# On Windows:
spacna_demo_env\Scripts\activate
# On macOS/Linux:
source spacna_demo_env/bin/activate

# Install dependencies
pip install torch==2.6.0+cu124 torchvision==0.21.0+cu124 --extra-index-url https://download.pytorch.org/whl/cu124
pip install opencv-python==4.12.0.88 numpy==2.2.6 pandas==2.3.3 scikit-learn==1.8.0 matplotlib==3.10.8 scipy==1.16.3
```


### 1.2 R Environment Setup

Open R or RStudio and run:

```r
# Install from CRAN
install.packages(c("Seurat", "parallelDist", "irlba", "ggplot2", "rootSolve", "patchwork", "glmnet", "dlm"))

# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install from Bioconductor
BiocManager::install(c("biomaRt", "ComplexHeatmap", "edgeR"))
```


## 2. Prepare Working Directory
Your directory structure should now look like:

```txt
demo_data
├── process_data
│   ├── gene_pos.rds
│   ├── seurat_object.rds
└── raw_data
    ├── spatial
    │   ├── aligned_fiducials.jpg
    │   ├── detected_tissue_image.jpg
    │   ├── exp_location.txt
    │   ├── spot.txt
    │   ├── tissue_hires_image.png
    │   ├── tissue_lowres_image.png
    │   ├── tissue_positions_list.csv
    │   └── scalefactors_json.json
    └── filtered_feature_bc_matrix.h5
```
## 3. Execute the Pipeline
### Step 1: Extract Image Features
```bash
# Make sure environment is activated
python get_spatial_feature.py --sample_dir "./demo_data/raw_data/spatial/" --radius 50
```
Verification: Check outputs were created:

```bash
ls ./demo_data/raw_data/spatial/resnet50/    # Should contain features_pca.txt
ls ./demo_data/raw_data/spatial/plot/        # Should contain contains graph images corresponding to different image thresholds. Based on the plot, select an appropriate image_thre for Step 2.
```

### Step 2: Detect CNAs with SpaCNA
```bash
Rscript SpaCNA.R \
    --sample_dir="./demo_data/process_data" \
    --image_dir="./demo_data/raw_data/spatial" \
    --plot_dir="./demo_data/result" \
    --normal_clusters="1,2" \
    --image_thre=0.2
```
**Note:** Adjust --normal_clusters based on your specific data. Check your seurat_object.rds for cluster information. Adjust --image_thre based on plots generated from Step 1.

**Verification:**
```bash
ls ./demo_data/results/     # Should contain CNA results and plots
```

Your results should now look like:

```txt
demo_data/
└──result/
│  ├── bk_bic.rds
│  ├── cns_clone_annot.png
│  ├── cns_init.png
│  ├── cns.png
│  ├── cns.rds
│  ├── copy_ratio.png
│  ├── count_norm.rds
│  └── hrmf_para.rds
```

### Step 3: Estimate Tumor Content (Optional)
```bash
Rscript estimate_tumor_content.R \
    --sample_dir="./demo_data/process_data" \
    --plot_dir="./demo_data/result" \
    --spacna_dir="./demo_data/result" \
    --K=7 \
    --cnv_thre==10
```
**Verification:**
```bash
ls ./demo_data/results/     # Should contain CNA results and plots
```
Your results should contain: spacna_tumor_content.png and tumor_content_df.rds

### Step 4: Detect Tumor Edge (Optional)

```bash
Rscript estimate_tumor_edge.R \
--sample_dir="./demo_data/process_data" \
    --plot_dir="./demo_data/result" \
    --spacna_dir="./demo_data/result" \
    --tumor_content_dir="./demo_data/result" \   
    --edge_thre=0.03 \
    --tumor_thre = 0.55
```

**Verification:**
```bash
ls ./demo_data/results/     # Should contain CNA results and plots
```
Your results should contain: spacna_edge.png and spacna_edge_score.png

