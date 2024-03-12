# Need to load the trace first
library(readr)
library(tidyverse)
library(lubridate)

# FUnctions to use
toSeconds <- function(x) {
  if(grepl("m", x)) {
    parts <- unlist(strsplit(x, "m|s"))
    return(as.numeric(parts[1]) * 60 + as.numeric(parts[2]))
  } else {
    return(as.numeric(gsub("s", "", x)))
  }
}

# Params
trace_txt <- here::here("run_poster/trace-20240311-71674714.txt")
# The data
raw_df <- read_tsv(trace_txt)
df <- raw_df %>% 
      mutate(process = tolower(process)) %>%
      select(-memory) %>%
  filter(str_detect(process, "cv")) %>%
  mutate(process = str_extract(process, "(mogonet|diablo|cooperative_learning)_(preprocess|train|predict)")) %>%
mutate(process = str_replace(process, "cooperative_learning", "cooperative-learning")) %>%
  separate(process, into = c("method", "action"), sep = "_", extra = "merge") %>%
  mutate(raw_seconds = sapply(duration, toSeconds) |> as.numeric()) %>%
  mutate(method = str_replace(method, "cooperative-learning", "cooperative_learning")) %>%
  filter(!is.na(method))


df %>%
  write.csv("transformed_time.csv", row.names = FALSE)

# Colors
custom_colors <- c("diablo" = '#1f77b4', 
                   "cooperative_learning" = '#ff7f0e', 
                   "mogonet" =  '#2ca02c')


df %>%
  ggplot(aes(x=method, y=raw_seconds, fill=method)) +
  stat_boxplot(geom ='errorbar', width=0.5) +
  geom_boxplot(outlier.color = "red") + 
  scale_fill_manual(values=custom_colors) +
  scale_y_log10() +
  labs(x = "Method", y = "Log10 of computation time in seconds") +
  theme_bw()


