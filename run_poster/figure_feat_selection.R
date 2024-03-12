# Use this to visualize feature selcetion stuff
library(here)
library(tidyverse)
library(UpSetR)
library(ggupset)
# Select feat
selectFeat <- function(df, method) {
  df %>%
    dplyr::filter(method == {{method}}) %>%
    dplyr::pull(feat)
}
# Combining many csvs
combine_csvs <- function(p) {
  all_files <- list.files(here(p))
  # allocate mem
  to_bind <- list()
  for (csv in all_files) {
    full_path <- here(p, csv)
    df <- read.csv(full_path)
    to_bind[[csv]] <- df
  }
  all_df <- bind_rows(to_bind) %>% as_tibble()
  return(all_df)
}
# Parameters used
mogonet_p <- here("features_selection_results/mogonet/")
diab_p <- here("features_selection_results/diablo/")
cplr_p <- here("features_selection_results/cplr/")


# Load those dfs
diab_dfs <- combine_csvs(diab_p)
cplr_dfs <- combine_csvs(cplr_p) %>% rename(feat=feature)
mogonet_dfs <- combine_csvs(mogonet_p) %>%
               mutate(dataset_name = ifelse(dataset_name == "multiomics-sim3", 
                                            "GSE123", dataset_name)
                      ) %>%
               select(-imp) %>%
               rename(feat = feat_name, view = omics) %>%
              # Massive check to rename
               select(view, feat, method, dataset_name) %>%
              mutate(view = case_when(
                dataset_name %in% c("GSE123", "multiomics-sim1", "multiomics-sim2") ~
                  ifelse(view == 0, "B_BLOCK", 
                         ifelse(view == 1, "A_BLOCK", "C_BLOCK")
                         ),
                dataset_name == "rosmap" ~
                  ifelse(grepl("ENSG", feat), "genomics", 
                         ifelse(grepl("cg", feat), "epigenomics", "transcriptomics")
                  ),
                dataset_name %in%  c("tcga", "GSE71669") ~
                  ifelse(grepl("hsa-", feat), "mirna", "unknown")
                  )
              )

# add another one to combine it for late usage
all_df <- bind_rows(cplr_dfs, diab_dfs, mogonet_dfs)

# Makes it to list
a_list <- list(
  COOPERATIVE_LEARNING = selectFeat(all_df, "cplr"),
  DIABLO = selectFeat(all_df, "diablo"),
  MOGONET = selectFeat(all_df, "mogonet")
               )

# reference https://dethwench.com/making-upset-plots-in-r-with-upsetr/

#I created different options so I could see which 
#ones I liked best
text_scale_options1 <- c(1, 1, 1, 1, 0.75, 1)
text_scale_options2 <- c(1.3, 1.3, 1, 1, 2, 0.75)
text_scale_options3 <- c(1.5, 1.25, 1.25, 1, 2, 1.5)

#setting colors
#this can also be done with hexadecimal
main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")
# Other stuff
#set_vars <- c("DIABLO", "MOGONET")
# PLot it
up_p <- UpSetR::upset(fromList(a_list),
      order.by = "freq",
      empty.intersections = "on",
      mainbar.y.label = "Total number of selected features",
      sets.x.label = "Number of overlapped features selected",
      show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      text.scale=text_scale_options3,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col
      )

#up_p

p_dat <- all_df %>% 
  distinct(method, feat, .keep_all=TRUE) %>%
  mutate(method = case_when(
    method == "cplr" ~ "COOPERATIVE_LEARNING",
    method == "mogonet" ~ "MOGONET",
    method == "diablo" ~ "DIABLO"
         )
    ) %>%
  group_by(feat) %>%
  summarize(method = list(method)) 

# plot params
width <- 12
height <- 12
fontsize <- 12

upset_plot <- p_dat %>%
  ggplot(aes(x=method)) +
  geom_bar(width = 0.5, fill = "skyblue") +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(order_by = "degree") +
  theme_bw() +
  theme(axis.title.y = element_text(vjust=-50)) +
  theme_combmatrix(combmatrix.label.text = element_text(size=fontsize),
                  combmatrix.panel.point.size = 5,
                  combmatrix.label.extra_spacing = fontsize
                  ) +

  labs(x="Method", y = "Numbers of chosen features")

upset_plot




ggsave("upset_plot.png", upset_plot, width=width, height=height)




