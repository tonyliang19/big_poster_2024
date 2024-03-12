library(here)
library(MultiAssayExperiment)
library(mixOmics)
library(multiview)
library(tidyr)
library(dplyr)
library(magrittr)
tcga_path <- here("data/real_data/breast_tcga/")
rosmap_path <- here("data/real_data/rosmap/")


# To binary factor
toBinFct <- function(x) {
  # Try to see if its numeric or string
  if (!is.numeric(x)) {
    x <- ifelse(x == "yes", 1, 0)
  }
  x <- factor(x, levels = c(0, 1))
  return(x)
}

# Read the data in first 
tcga_mae <- loadHDF5MultiAssayExperiment(paste0(tcga_path, "breast_tcga_mae_data"))
rosmap_mae <- loadHDF5MultiAssayExperiment(paste0(rosmap_path, "mae_data"))
# Split them to X and y
tcga_X <- tcga_mae@ExperimentList@listData |> lapply(t) |> lapply(as.matrix)
tcga_y <- tcga_mae$response |> toBinFct()
rosmap_X <- rosmap_mae@ExperimentList@listData |> lapply(t) |> lapply(as.matrix)
rosmap_y <- rosmap_mae$response |> toBinFct()

# So they have columns as the patients/samples names, already transposed
# so thats why checking rownames
lapply(tcga_X, rownames)
lapply(rosmap_X, rownames)

# Then you want to run over the model (one shot)
# TCGA
createKeepX <- function(x_list, ncomp=2) {
  # For each omics, only get 10% of total features
  # Assume the input list is already after transforming to rows being common
  # samples
  keep_list <- lapply(x_list, function(x) {
    feat_num <- floor(ncol(x) / 10)
    # repat the number of features for ncomp times as a vector
    feat_vec <- rep(feat_num, ncomp)
    return(feat_vec)
  })
  # assign to have the matching omics names
  names(keep_list) <- names(x_list)
  return(keep_list)
}

# Should fit model

extract_names <- function(var_list) {
  selected_feat_list <- lapply(var_list, function(omic_list) {
    # just retrive the name column and not the value col
    return(omic_list$name)
  })
  return(selected_feat_list)
}



diablo_tcga <- block.splsda(tcga_X, tcga_y, keepX=createKeepX(tcga_X))
cplr_tcga <- multiview(tcga_X, tcga_y, family = binomial(), rho = 0)

# Rosmap
diablo_rosmap <- block.splsda(rosmap_X, rosmap_y, keepX=createKeepX(rosmap_X))
cplr_rosmap <- multiview(rosmap_X, rosmap_y, family = binomial(), rho = 0)
cved <- cv.multiview(rosmap_X, rosmap_y, family = binomial(), type.measure = "deviance", rho=0.5)


# diablo get features
tcga_diab_var_list <- selectVar(diablo_tcga)
rosmap_diab_var_list <- selectVar(diablo_rosmap)
# Last column is ncomp, so discard that part, hence -1
tcga_diab_feats <- extract_names(tcga_diab_var_list[1:length(tcga_diab_var_list) - 1])
rosmap_diab_feats <- extract_names(rosmap_diab_var_list[1:length(rosmap_diab_var_list) - 1])


purrr::map_df(tcga_diab_feats, ~ tibble(col1 = names(tcga_diab_feats), col2 = .x))


aaa <- tibble::enframe(tcga_diab_feats, name="view") %>%
  unnest(value) %>%
  rename(feat=value)




t# cplr get features
# So from the coefficients get those non 0
getFeatCPLR <- function(fit) {
  # Needs to be multiview fit
  s <- fit$lambda.min
  feats_df <- coef_ordered(fit, s=s) %>%
    as_tibble() %>%
    select(view, view_col) %>%
    arrange(view) %>%
    rename(feature=view_col)
  return(feats_df)
}

# This how it calls cplr instead
cved <- cv.multiview(rosmap_X, rosmap_y, family = binomial())
s = cved$lambda.min
jjj <- coef_ordered(cved, s=s) %>%
  as_tibble() %>%
  select(view, view_col) %>%
  rename(feature=view_col)

# Main funs to implement

main_diablo <- function(path) {
  # Read the data in
  # fit the model
  
  # and need to keepX (requires rows are samples and not cols)
  
  # then from its selectedVars extract features as list
  # lastly convert the list into a dataframe of
  # view (omics), feature
  
}


main_cplr <- function(path) {
  # Read the data in
  # fit the model (but with cv) ....
  # TODO: need to fix this above
  
  # Then get an s likely to be lambda.min or lambda1se if preffered
  # but 1se provides very little feats
  # TODO: add stuff
  
  # then could use coef_ordered to select relevant cols
  # since it would include the importance val as well
  # lastly sort it by view (omics)
  # and rename view_col to be feature
  # Already a dataframe so good!
  
}

