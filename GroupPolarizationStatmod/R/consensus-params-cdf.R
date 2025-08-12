library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read data
df <- read_csv("data/StudiesAnalysis.csv")

# Filter to plausibly spurious & included
df_plot <- df %>%
  filter(PlausibleInt == 1, IncludeInt == 1) %>%
  select(ArticleTag, LatentSDPre, LatentSDPost) %>%
  pivot_longer(cols = c(LatentSDPre, LatentSDPost),
               names_to = "phase", values_to = "latent_sd") %>%
  mutate(phase = ifelse(phase == "LatentSDPre", "Pre", "Post"))

# Base plot object
p <- ggplot() +
  # Per-article ECDFs (thin + transparent)
  stat_ecdf(data = df_plot, aes(x = latent_sd, group = interaction(ArticleTag, phase),
                                color = phase, linetype = phase),
            geom = "step", alpha = 0.2, size = 0.4) +
  # Overall ECDFs (bold)
  stat_ecdf(data = df_plot, aes(x = latent_sd, color = phase, linetype = phase),
            geom = "step", size = 1.2) +
  scale_color_manual(values = c("Pre" = "#1b9e77", "Post" = "#d95f02")) +
  labs(x = "Latent standard deviation",
       y = "Cumulative fraction of conditions",
       title = "Pre vs. Post deliberation latent SDs (plausibly spurious cases)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

p


library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read data
df <- read_csv("data/StudiesAnalysis.csv")

# Filter to plausibly spurious and included
df_plot <- df %>%
  filter(PlausibleInt == 1, IncludeInt == 1) %>%
  select(TreatmentTag, ArticleTag, LatentSDPre, LatentSDPost) %>%
  pivot_longer(cols = c(LatentSDPre, LatentSDPost),
               names_to = "Phase", values_to = "latent_sd") %>%
  mutate(Phase = ifelse(Phase == "LatentSDPre", "Pre", "Post"))

# Plot CDF
p <- ggplot(df_plot, aes(x = latent_sd, color = Phase, linetype = Phase)) +
  stat_ecdf(size = 1) +
  labs(x = "Latent standard deviation",
       y = "Cumulative fraction of conditions") + #,
       #title = "Pre vs. Post deliberation latent SDs (plausibly false cases)") +
  scale_color_manual(values = c("Pre" = "#1b9e77", "Post" = "#d95f02")) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.5, 0.8, 1.0),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 14) + 
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "lightgrey") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "lightgrey") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "lightgrey") +
  theme(panel.grid.minor = element_blank()) 
  
ggsave("paper/Figures/3.1/latent-sd-cdf.png", p, width = 6, height = 3.5, units = "in", dpi = 300)
  