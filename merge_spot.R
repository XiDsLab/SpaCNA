suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

#' Merge spatial spots by integer grid scaling
#'
#' @param seu_obj Seurat object with spatial coordinates
#' @param min_spots Minimum spots after merging (default: 1000)
#' @param max_spots Maximum spots after merging (default: 4000)
#' @param verbose Print progress messages (default: TRUE)
#' @return Merged Seurat object with grid coordinates
merge_spatial_spots <- function(seu_obj, min_spots = 1000, max_spots = 4000, verbose = TRUE) {
  
  # Get original data
  coords <- GetTissueCoordinates(seu_obj)
  expr_matrix <- GetAssayData(seu_obj)
  n_original <- nrow(coords)

  if (verbose) message(sprintf("Original spots: %d", n_original))
  
  # Check if merging is needed
  if (n_original <= max_spots) {
    if (verbose) message("No merging needed")
    coords_result <- get_coords_merge_obj(seu_obj, seu_obj)
    return(list(
      merged_obj = seu_obj,
      exp_location = coords_result$exp_location,
      spot_name = coords_result$spot_name
    ))
  }
  
  # Find best integer scaling factor
  scaling <- find_best_scaling(n_original, min_spots, max_spots)
  
  # Create new grid coordinates
  grid_data <- create_integer_grid(coords, scaling$k, verbose)
  
  # Assign spots to new grid cells
  assignments <- assign_to_integer_grid(coords, grid_data, scaling$k, verbose)
  
  # Merge expression matrix
  merged_expr <- merge_expression_integer(expr_matrix, assignments, scaling$target_spots, verbose)
  
  # Create new Seurat object
  merged_obj <- create_merged_seurat_obj(merged_expr, assignments$centers, assignments)

  # Get exp_location.txt and spot.txt
  coords_result <- get_coords_merge_obj(merged_obj,seu_obj)
  
  if (verbose) {
    final_spots <- ncol(merged_expr)
    message(sprintf("\nMerge complete: %d → %d spots", n_original, final_spots))
    message(sprintf("Grid: %dx%d cells", grid_data$grid_dim, grid_data$grid_dim))
  }
  
  return(list(
    merged_obj = merged_obj,
    exp_location = coords_result$exp_location,
    spot_name = coords_result$spot_name
  ))
}

#' Find best integer scaling factor (k^2) where k is integer
find_best_scaling <- function(n_original, min_spots, max_spots) {
  if (n_original <= max_spots) {
    return(list(k = 1, factor = 1, target_spots = n_original))
  }
  
  k <- ceiling(sqrt(n_original / max_spots))
  k <- max(k, 1)
  
  while (TRUE) {
    target_spots <- ceiling(n_original / k^2)
    
    if (target_spots > max_spots) {
      k <- k + 1 
    } else if (target_spots < min_spots && k > 1) {
      k <- k - 1
    } else {
      break
    }
  }
  
  return(list(
    k = k,
    factor = k^2,
    target_spots = target_spots
  ))
}

#' Create integer grid based on scaling factor
create_integer_grid <- function(coords, k, verbose) {
  x_range <- range(coords$x)
  y_range <- range(coords$y)
  
  estimated_target <- ceiling(nrow(coords) / k^2)
  grid_dim <- ceiling(sqrt(estimated_target))
  
  x_centers <- seq(x_range[1], x_range[2], length.out = grid_dim)
  y_centers <- seq(y_range[1], y_range[2], length.out = grid_dim)
  
  centers <- expand.grid(x = x_centers, y = y_centers)
  rownames(centers) <- sprintf("center_%d", 1:nrow(centers))
  
  if (verbose) {
    message(sprintf("Created %d×%d = %d merging centers", 
                   grid_dim, grid_dim, nrow(centers)))
  }
  
  return(list(
    centers = centers,
    grid_dim = grid_dim,
    k = k
  ))
}

#' Assign spots to centers within square region
assign_to_integer_grid <- function(coords, grid_data, k, verbose) {
  centers <- grid_data$centers
  original_centers_count <- nrow(centers)
  
  if (nrow(centers) > 1) {
    unique_x <- unique(centers$x)
    unique_y <- unique(centers$y)
    
    x_interval <- min(diff(sort(unique_x))) 
    y_interval <- min(diff(sort(unique_y))) 
  } else {
    x_range <- range(coords$x)
    y_range <- range(coords$y)
    x_interval <- diff(x_range)
    y_interval <- diff(y_range)
  }
  
  x_half_side <- x_interval / 2
  y_half_side <- y_interval / 2
  n_centers <- nrow(centers)
  
  x_min <- centers$x - x_half_side
  x_max <- centers$x + x_half_side
  y_min <- centers$y - y_half_side
  y_max <- centers$y + y_half_side

  assignments <- integer(nrow(coords))
  
  for (i in 1:nrow(coords)) {
    in_region <- which(
      coords$x[i] >= x_min & 
      coords$x[i] < x_max & 
      coords$y[i] >= y_min & 
      coords$y[i] < y_max
    )  
    
    if (length(in_region) == 1) {
      assignments[i] <- in_region
    } else if (length(in_region) == 0) {
      dists <- sqrt((centers$x - coords$x[i])^2 + (centers$y - coords$y[i])^2)
      assignments[i] <- which.min(dists)
    } else {
      dists <- sqrt((centers$x[in_region] - coords$x[i])^2 + 
                    (centers$y[in_region] - coords$y[i])^2)
      assignments[i] <- in_region[which.min(dists)]
    }
  }
  
  spots_per_center <- tabulate(assignments, nbins = n_centers)
  

  non_empty <- which(spots_per_center > 0)
  n_kept <- length(non_empty)
  
  if (n_kept == 0) {
    stop("No centers received spots. Check your data.")
  }
  

  old_to_new <- integer(n_centers)
  old_to_new[non_empty] <- 1:n_kept
  

  new_assignments <- integer(nrow(coords))
  valid_spots <- which(assignments %in% non_empty)
  new_assignments[valid_spots] <- old_to_new[assignments[valid_spots]]
  

  new_centers <- centers[non_empty, , drop = FALSE]
  new_spots_per_center <- spots_per_center[non_empty]
  
  rownames(new_centers) <- sprintf("center_%d", 1:n_kept)
  
  if (verbose) {
    n_empty <- n_centers - n_kept
    message(sprintf("  Original centers: %d", n_centers))
    message(sprintf("  Non-empty centers: %d", n_kept))
    message(sprintf("  Empty centers removed: %d", n_empty))
    message(sprintf("  Average spots per center: %.1f", mean(new_spots_per_center)))
  }
  
  return(list(
    assignments = new_assignments,
    spots_per_center = new_spots_per_center,
    centers = new_centers,
    valid_spots = valid_spots,
    original_centers_count = original_centers_count
  ))
}


#' Merge expression matrix
merge_expression_integer <- function(expr_matrix, assignment_result, n_target, verbose) {
  assignments <- assignment_result$assignments
  spots_per_center <- assignment_result$spots_per_center
  valid_spots <- assignment_result$valid_spots
  
  n_centers <- length(spots_per_center)
  n_genes <- nrow(expr_matrix)

  valid_expr <- expr_matrix[, valid_spots, drop = FALSE]
  valid_assignments <- assignments[valid_spots]
  
  merged_expr <- Matrix(0, nrow = n_genes, ncol = n_centers, sparse = TRUE)
  
  if (verbose) {
    message(sprintf("Merging %d spots into %d centers...", length(valid_spots), n_centers))
    #pb <- txtProgressBar(min = 0, max = n_centers, style = 3)
  }
  
  for (center in 1:n_centers) {
    #if (verbose && center %% 100 == 0) {
    #  setTxtProgressBar(pb, center)
    #}
    
    spot_indices <- which(valid_assignments == center)
    
    if (length(spot_indices) > 0) {
      if (length(spot_indices) == 1) {
        merged_expr[, center] <- valid_expr[, spot_indices]
      } else {
        merged_expr[, center] <- Matrix::rowSums(valid_expr[, spot_indices, drop = FALSE])
      }
    }
  }
  
  #if (verbose) close(pb)

  rownames(merged_expr) <- rownames(expr_matrix)
  colnames(merged_expr) <- sprintf("merged_%04d", 1:n_centers)
  
  if (verbose) {
    message(sprintf("\nCreated %d merged spots", n_centers))
  }
  
  return(merged_expr)
}


#' Create merged Seurat object
create_merged_seurat_obj <- function(merged_expr, centers, assignment_result) {
  
  # Use filtered centers
  new_coords <- centers 
  rownames(new_coords) <- colnames(merged_expr)
  
  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = merged_expr,
    project = 'SlideSeq',
    assay = "Spatial"
  )
  
  spatial_coords <- new_coords[, c("x", "y")]
  colnames(spatial_coords) <- c("imagerow", "imagecol")
  obj[['image']] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = spatial_coords
  )

  colnames(obj) <- gsub("_","-",colnames(obj))
  
  return(obj)
}

get_coords_merge_obj <- function(merged_obj,seu_obj) {
  colnames(merged_obj) <- gsub("_","-",colnames(merged_obj))
  scale_factor <- seu_obj@images[[1]]@scale.factors$hires
  coordinates <- GetTissueCoordinates(merged_obj)
  pixel_coords <- coordinates[, c("x", "y")] * scale_factor
  spot_name <- rownames(coordinates)
   
  return(list(
    exp_location = pixel_coords,
    spot_name = spot_name
  ))
}


#' Command Line Interface for merge spatial spots
#' 
#' Usage: Rscript combine_spot.r --input <input.rds> --output <output_dir> [options]
#' 
#' @param --input       Path to input Seurat RDS file
#' @param --output      Path to output directory
#' @param --min-spots   Minimum spots after merging (default: 1000)
#' @param --max-spots   Maximum spots after merging (default: 4000)
#' @param --help        Show this help message
run_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    cat("
Usage: Rscript combine_spot.r --input <input.rds> --output <output_dir> [options]

Required:
  --input FILE       Input Seurat RDS file with spatial data
  --output DIR       Output directory (merged.rds and txt files will be saved here)

Options:
  --min-spots N      Minimum spots after merging (default: 1000)
  --max-spots N      Maximum spots after merging (default: 4000)
  -h, --help         Show this help message

Examples:
  Rscript combine_spot.r --input data/seurat_obj.rds --output ./output
  Rscript combine_spot.r --input data/seurat_obj.rds --output ./output --min-spots 1500 --max-spots 6000
")
    return(invisible())
  }
  
  get_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) > 0 && idx < length(args)) {
      return(args[idx + 1])
    }
    return(NULL)
  }
  
  input_file <- get_arg("--input")
  output_dir <- get_arg("--output")
  
  min_spots_arg <- get_arg("--min-spots")
  min_spots <- if (is.null(min_spots_arg)) 1000 else as.integer(min_spots_arg)
  
  max_spots_arg <- get_arg("--max-spots")
  max_spots <- if (is.null(max_spots_arg)) 4000 else as.integer(max_spots_arg)
  
  verbose <- TRUE
  
  if (is.null(input_file)) {
    stop("Error: --input is required. Use --help for usage information.")
  }
  
  if (is.null(output_dir)) {
    stop("Error: --output is required. Use --help for usage information.")
  }
  
  if (!file.exists(input_file)) {
    stop(sprintf("Error: Input file not found: %s", input_file))
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Created output directory: %s", output_dir))
  }
  
  base_name <- gsub("\\.rds$", "", basename(input_file))
  output_rds <- file.path(output_dir, paste0(base_name, "_merged.rds"))
  
  message("Loading Seurat object...")
  seu_obj <- readRDS(input_file)
  
  message(sprintf("  Loaded: %d spots, %d genes", ncol(seu_obj), nrow(seu_obj)))
  
  message("\nMerging spots...")
  start_time <- Sys.time()
  
  result <- merge_spatial_spots(
    seu_obj = seu_obj,
    min_spots = min_spots,
    max_spots = max_spots,
    verbose = verbose 
  )

  merged_obj <- result$merged_obj
  exp_location <- result$exp_location
  spot_name <- result$spot_name
  
  message("\nSaving results...")
  
  saveRDS(merged_obj, file = output_rds)
  message(sprintf("  Saved: %s", output_rds))
  
  if (!is.null(exp_location) && !is.null(spot_name)) {
    exp_location_file <- file.path(output_dir, "exp_location.txt")
    spot_name_file <- file.path(output_dir, "spot.txt")
    
    write.table(exp_location, 
              file = exp_location_file,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
    message(sprintf("  Saved: %s", exp_location_file))

    write.table(spot_name,
              file = spot_name_file,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
    message(sprintf("  Saved: %s", spot_name_file))
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  cat("\n========================================\n")
  cat("Complete\n")
  cat("========================================\n")
  cat(sprintf("Original spots:  %d\n", ncol(seu_obj)))
  cat(sprintf("Merged spots:    %d\n", ncol(merged_obj)))
  cat(sprintf("Elapsed time:    %.1f seconds\n", elapsed))
  cat(sprintf("Output directory: %s\n", output_dir))
  cat("========================================\n")
}

if (sys.nframe() == 0 && !interactive()) {
  run_cli()
}