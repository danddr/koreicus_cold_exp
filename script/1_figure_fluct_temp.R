setwd("/home/dared/GitHub/koreicus_cold_exp/")
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
#set out directory
outdir<-"outputs/"

df <- read_xlsx("data/temp_setting.xlsx")

# Fill missing values with interpolation (to keep lines continuous)
df_interp <- df %>%
  mutate(across(starts_with("5_to"), ~approx(Hours, ., Hours, rule = 2)$y))

# Reshape for ggplot
df_long <- df_interp %>%
  pivot_longer(cols = -Hours, names_to = "Treatment", values_to = "Temperature")

# Plot
p <- ggplot(df_long, aes(x = Hours, y = Temperature, color = Treatment)) +
  geom_line(size = 1.2) +
  # geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype="dashed")+
  # scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "black")) +
  scale_color_manual(
    values = c(
      "5_to_0"   = "#bdd7e7",  # darker light blue
      "5_to_-5"  = "#6baed6",  # medium-light blue
      "5_to_-10" = "#4292c6",  # medium blue
      "5_to_-15" = "#2171b5",  # dark blue
      "5_to_-20" = "#084594",  # very dark blue
      "Control"  = "#d95f02"   # contrast orange  
      
    ),
    labels = c(
      "5_to_0"   = "5°C → 0°C",
      "5_to_-5"  = "5°C → −5°C",
      "5_to_-10" = "5°C → −10°C",
      "5_to_-15" = "5°C → −15°C",
      "5_to_-20" = "5°C → −20°C",
      "Control"  = "Control (5°C constant)"
    ),
    breaks = c(
      "5_to_0",
      "5_to_-5",
      "5_to_-10",
      "5_to_-15",
      "5_to_-20",
      "Control"
    ))+
  ylim(-25, 10)+
  scale_x_continuous(breaks = seq(0, 24, 4)) +
    labs(
      x = "Time (hours)",
    y = "Temperature (°C)",
    color = NULL
  ) +
  coord_fixed(ratio = 1)+
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x  = element_text(size = 30),
    axis.text.y  = element_text(size = 30),
  )
p
outname <- paste0(outdir, "/fluct_temp",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 15, height = 10, device='png', dpi=320)


