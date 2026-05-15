# Setting domestic CRAN and Bioconductor mirrors ----
options("repos" = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor/")

# Installing required packages ----
if(!"devtools" %in% installed.packages()){install.packages("devtools")}
if(!"httr" %in% installed.packages()){install.packages("httr")}
if(!"bio3d" %in% installed.packages()){install.packages("bio3d")}
if(!"openxlsx" %in% installed.packages()){install.packages("openxlsx")}
if(!"dplyr" %in% installed.packages()){install.packages("dplyr")}
if(!"gluedocking" %in% installed.packages()){devtools::install_github("RightSZ/gluedocking")}

# Loading required packages ----
options(warn = -1)
library(gluedocking)
library(httr)
library(bio3d)
library(openxlsx)
library(dplyr)

# Step 1: Prepare configuration files ----
prepare_for_gluedocking(
  python_path = "D:\\Docks\\autodockTools\\python.exe",
  prepare_receptor_script = "D:\\Docks\\autodockTools\\Lib\\site-packages\\AutoDockTools\\Utilities24\\prepare_receptor4.py",
  prepare_ligand_script = "D:\\Docks\\autodockTools\\Lib\\site-packages\\AutoDockTools\\Utilities24\\prepare_ligand4.py",
  prepare_split_alt_script = "D:\\Docks\\autodockTools\\Lib\\site-packages\\AutoDockTools\\Utilities24\\prepare_pdb_split_alt_confs.py",
  obabel_path = "D:\\Docks\\OpenBabel\\OpenBabel-3.1.1\\obabel.exe",
  vina_path = "D:\\Docks\\Miniconda3\\Library\\bin\\qvinaw.exe",
  force = T
)

# Step 2: Download PDB files of protein receptors ----
if(!dir.exists('result/receptors')){dir.create('result/receptors', recursive = T)}
pdb_files <- lapply(list.files('result/receptors'), function(f){
  paste0('result/receptors/',f)
})

# Step 3: Download SDF files of compound ligands ----
if(!dir.exists('result/ligands')){dir.create('result/ligands', recursive = T)}
ligand_files <- lapply(list.files('result/ligands'), function(f){
  paste0('result/ligands/',f)
})

# Step 4: Preprocess protein receptor data (removing water, hydrogen, original ligands, etc.) ----
# Check if preprocessed files already exist
trimmed_dir <- "result/trimmed"
if (dir.exists(trimmed_dir) && length(list.files(trimmed_dir, pattern = "\\.pdb$")) > 0) {
  message("Trimmed files already exist. Importing directly...")
  trimmed_files <- list.files(trimmed_dir, pattern = "\\.pdb$", full.names = TRUE)
} else {
  trimmed_files <- trim_receptor("result/receptors", output_dir = trimmed_dir)
}

# Check split_alt files
split_dir <- "result/split_alt"
if (dir.exists(split_dir) && length(list.files(split_dir, pattern = "\\.pdb$")) > 0) {
  message("Split_alt files already exist. Importing directly...")
  pbd_files <- list.files(split_dir, pattern = "\\.pdb$", full.names = TRUE)
} else {
  pbd_files <- split_alt(
    inputs = trimmed_files,
    output_dir = split_dir,
    keep_label = "A"
  )
}

# 4.1 Convert PDB to PDBQT format ----
receptor_dir <- "result/prepared_receptors"
if (dir.exists(receptor_dir) && length(list.files(receptor_dir, pattern = "\\.pdbqt$")) > 0) {
  message("Receptor PDBQT files already exist. Importing directly...")
  receptor_pdbqt <- list.files(receptor_dir, pattern = "\\.pdbqt$", full.names = TRUE)
} else {
  receptor_pdbqt <- prepare_receptor(
    pbd_files,
    output_dir = receptor_dir
  )
}

# Step 5: Preprocess compound ligand data ----
# 5.1 Convert SDF to MOL2 format ----
mol2_dir <- "result/ligands_mol2"
if (dir.exists(mol2_dir) && length(list.files(mol2_dir, pattern = "\\.mol2$")) > 0) {
  message("MOL2 files already exist. Importing directly...")
  converted_files <- list.files(mol2_dir, pattern = "\\.mol2$", full.names = TRUE)
} else {
  converted_files <- convert_molecule(
    input_file = "result/ligands",
    output_format = "mol2",
    output_dir = mol2_dir
  )
}

# 5.2 Convert to PDBQT format ----
ligand_dir <- "result/prepared_ligands"
if (dir.exists(ligand_dir) && length(list.files(ligand_dir, pattern = "\\.pdbqt$")) > 0) {
  message("Ligand PDBQT files already exist. Importing directly...")
  ligand_pdbqt <- list.files(ligand_dir, pattern = "\\.pdbqt$", full.names = TRUE)
} else {
  ligand_pdbqt <- prepare_ligand(
    converted_files,
    output_dir = ligand_dir
  )
}

# Step 6: Calculate active docking pocket region ----
# Check whether recalculation of box parameters is needed (adjustable as needed)
box_file <- "result/box_params.rds"
if (file.exists(box_file)) {
  message("Box parameter file already exists. Importing directly...")
  box_params <- readRDS(box_file)
} else {
  box_params <- calculate_box(receptor_pdbqt, padding = 5)
  # Save box parameters for subsequent use
  saveRDS(box_params, box_file)
}
print(box_params)

# Step 7: Generate configuration files for molecular docking ----
config_dir <- "result/configs"

# Preliminary data size check
cat("\n========== Data Size Check ==========\n")
cat("Number of receptor PDBQT files:", length(receptor_pdbqt), "\n")
cat("Number of ligand PDBQT files:", length(ligand_pdbqt), "\n")
cat("Number of box parameter rows:", nrow(box_params), "\n")
expected_count <- length(receptor_pdbqt) * length(ligand_pdbqt)
cat("Expected number of configuration files:", expected_count, "\n")
cat("======================================\n")

# Check existing configuration files
if (!dir.exists(config_dir)) {
  dir.create(config_dir, recursive = TRUE)
  existing_files <- character(0)
  message("Created configuration folder: ", config_dir)
} else {
  existing_files <- list.files(config_dir, pattern = "\\.txt$", full.names = TRUE)
}

# Analyze generated and missing files
if (length(existing_files) > 0) {
  cat("\n========== Configuration File Check ==========\n")
  cat("Existing configuration files:", length(existing_files), "\n")
  
  # Extract receptor and ligand information from existing filenames to identify which combinations have been generated
  # Assuming filename format: receptor_ligand_config.txt
  existing_pairs <- data.frame(
    file = existing_files,
    receptor = gsub("_(.*?)_config\\.txt$", "", basename(existing_files)),
    ligand = gsub("^.*?_(.*?)_config\\.txt$", "\\1", basename(existing_files))
  )
  
  # Generate all possible combinations
  all_combinations <- expand.grid(
    receptor = basename(receptor_pdbqt),
    ligand = basename(ligand_pdbqt),
    stringsAsFactors = FALSE
  )
  all_combinations$file_name <- paste0(all_combinations$receptor, "_", 
                                       all_combinations$ligand, "_config.txt")
  
  # Identify missing combinations
  generated_files <- basename(existing_files)
  missing_combinations <- all_combinations[!all_combinations$file_name %in% generated_files, ]
  
  cat("Number of generated combinations:", nrow(all_combinations) - nrow(missing_combinations), "\n")
  cat("Number of missing combinations:", nrow(missing_combinations), "\n")
  
  if (nrow(missing_combinations) > 0) {
    cat("\nExamples of missing receptor-ligand combinations:\n")
    print(head(missing_combinations[, c("receptor", "ligand")], 5))
    
    # Prompt whether to generate missing files
    response <- readline(prompt = "\nGenerate missing configuration files? (y/n): ")
    
    if (tolower(response) == "y") {
      cat("\nInitiating generation of missing configuration files...\n")
      
      # Process missing combinations in batches
      missing_receptors <- unique(missing_combinations$receptor)
      missing_ligands <- unique(missing_combinations$ligand)
      
      # Call write_configs function only for missing combinations
      # Note: Adjustments based on actual conditions are required here, as write_configs may need complete receptor and ligand paths
      # Approach 1: If write_configs supports specifying specific combinations
      config_files <- write_configs(
        receptor_paths = receptor_pdbqt[basename(receptor_pdbqt) %in% missing_receptors],
        ligand_paths = ligand_pdbqt[basename(ligand_pdbqt) %in% missing_ligands],
        box_df = box_params,
        output_dir = config_dir,
        exhaustiveness = 1
      )
      
      cat("Generation of missing configuration files completed. Number of new files added:", length(config_files), "\n")
    } else {
      cat("Skipped generation of missing configuration files\n")
    }
  } else {
    cat("\nAll configuration files have been completely generated. No further action required!\n")
  }
  
  # Final statistics
  final_files <- list.files(config_dir, pattern = "\\.txt$", full.names = TRUE)
  cat("\n========== Final Status ==========\n")
  cat("Total number of configuration files:", length(final_files), "\n")
  cat("Expected number:", expected_count, "\n")
  
  if (length(final_files) == expected_count) {
    cat("✓ Configuration files have been completely generated!\n")
  } else {
    cat("⚠ Configuration files are incomplete. Missing", expected_count - length(final_files), "files\n")
  }
  cat("===================================\n")
  
  config_files <- final_files
  
} else {
  # No existing files, generate all from scratch
  cat("\nNo existing configuration files found. Generating all files...\n")
  cat("Approximately", expected_count, "configuration files will be generated. This may take some time...\n")
  
  response <- readline(prompt = "Proceed? (y/n): ")
  if (tolower(response) == "y") {
    config_files <- write_configs(
      receptor_paths = receptor_pdbqt,
      ligand_paths = ligand_pdbqt,
      box_df = box_params,
      output_dir = config_dir,
      exhaustiveness = 1
    )
    cat("Number of generated configuration files:", length(config_files), "\n")
  } else {
    stop("User cancelled configuration file generation")
  }
}

# Step 8: Execute AutoDock Vina for molecular docking ----
results <- run_vina(
  config_paths = "result/configs",
  output_dir = "result/docked",
  logs_dir = "result/logs",
  cpu = 8
)

# Step 9: Evaluate docking results; re-run failed docking tasks ----
check_logs(
  logs_dir = "result/logs",
  output_dir = "result/docked",
  config_paths = "result/configs"
)

# Step 10: Organize docking results ----
results_df <- parse_logs(logs_dir = "result/logs")

write.xlsx(results_df, "result/docking_results_all.xlsx", rowNames = F)





