## Installing ggplot2 library in R if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

## Here, we create a function to generate plots for each nucleotidic pair using ggplot
generate_plots <- function(pair_df, pair_name) {
  pair_df2 <- data.frame(
    Distance = seq(1, 20),
    Score = unlist(pair_df[, 3:22])
  )
  
  ggplot(pair_df2, aes(x = Distance, y = Score)) +
    geom_line() +
    labs(title = paste("Pair", pair_name),
         x = "Sequential Distance",
         y = "Pseudo-Energy Score") +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white"))
}

## Uploading the CSV files that we obtained from the obj_func_training.py
df <- read.csv("pseudo_energy_data.csv")

## Then we use unique to take into account only the unique nucleotide pairs
unique_pairs <- unique(df[, c("Nucleotide.1", "Nucleotide.2")])

## We then use a foor loop to generate plots for each pair
for (i in 1:nrow(unique_pairs)) {
  pair <- unique_pairs[i, ]
  pair_df <- subset(df, Nucleotide.1 == pair$Nucleotide.1 & Nucleotide.2 == pair$Nucleotide.2)
  plot <- generate_plots(pair_df, paste(pair$Nucleotide.1, pair$Nucleotide.2, sep = '-'))
  plot_name <- paste("Pair", pair$Nucleotide.1, pair$Nucleotide.2, "plot.png", sep = '_')
  ggsave(plot_name, plot)
}
