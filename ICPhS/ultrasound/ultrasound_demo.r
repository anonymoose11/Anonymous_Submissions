# Demonstration script
#The analysis code was written by the authors; generative AI was used only to standardize formatting and improve consistency in presentation.

# -------------------------------------------------------------------------
# 1. Libraries
# -------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)
library(mgcv)
library(itsadug)

# -------------------------------------------------------------------------
# 1. Load and normalize Spanish data
# -------------------------------------------------------------------------

dat <- read_excel("ultrasound_demo.xlsx") %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  ) %>%
  group_by(participant) %>%
  mutate(
    z_x = (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE),
    z_y = (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
  ) %>%
  ungroup()

# -------------------------------------------------------------------------
# 2. Remove contours containing outlying points
# -------------------------------------------------------------------------

ref <- dat %>%
  group_by(participant, vowel, closest_angle) %>%
  summarise(
    med_y = median(z_y, na.rm = TRUE),
    mad_y = mad(z_y, constant = 1.4826, na.rm = TRUE),
    .groups = "drop"
  )

bad_files <- dat %>%
  left_join(ref, by = c("participant", "vowel", "closest_angle")) %>%
  mutate(
    mad_y = ifelse(is.na(mad_y) | mad_y < 1e-6, NA_real_, mad_y),
    outlier = !is.na(mad_y) & abs(z_y - med_y) > 3 * mad_y
  ) %>%
  group_by(participant, vowel, file) %>%
  summarise(drop = any(outlier), .groups = "drop") %>%
  filter(drop)

sp <- dat %>%
  anti_join(bad_files, by = c("participant", "vowel", "file"))

write_xlsx(sp, "spanish_cleaned.xlsx")

# -------------------------------------------------------------------------
# 3. Descriptive contour plot
# -------------------------------------------------------------------------

vowel_colours <- c(
  "i" = "#A869C7", "e" = "#FA5252",
  "u" = "#2B67E0", "o" = "#4A9663", "a" = "#F76F00"
)

sp_plot <- ggplot(sp, aes(z_x, -z_y, colour = vowel)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.5) +
  scale_x_reverse() +
  scale_colour_manual(values = vowel_colours) +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(family = "sans", size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "← Front",
    y = "Height →",
    colour = "Vowel",
    title = "Spanish normalized tongue contours"
  )

ggsave(
  "spanish_normalized_contours.png",
  sp_plot, width = 8, height = 6, dpi = 300
)

# -------------------------------------------------------------------------
# 4. GAMMs
# -------------------------------------------------------------------------

prepare_pair <- function(data, vowels, reference) {
  data %>%
    filter(vowel %in% vowels) %>%
    mutate(
      z_y = -z_y,
      vowel = relevel(factor(vowel), reference),
      participant = factor(participant),
      file = factor(file),
      word = factor(word),
      closest_angle = as.numeric(closest_angle)
    )
}

fit_pair <- function(data) {

  # Initial model
  m1 <- bam(
    z_y ~
      s(participant, bs = "re") +
      s(file, bs = "re") +
      s(closest_angle, participant, bs = "fs"),
    data = data,
    method = "fREML",
    discrete = TRUE
  )
  rho1 <- acf_resid(m1)[2]

  # Define the start of each autocorrelation series
  data <- data %>%
    arrange(participant, word, closest_angle) %>%
    group_by(participant, word) %>%
    mutate(AR.start = closest_angle == min(closest_angle)) %>%
    ungroup()

  # Intermediate AR(1) model
  m2 <- bam(
    z_y ~
      s(closest_angle) +
      s(participant, bs = "re") +
      s(file, bs = "re") +
      s(closest_angle, participant, bs = "fs"),
    AR.start = data$AR.start,
    rho = rho1,
    data = data,
    method = "fREML",
    discrete = TRUE
  )
  rho2 <- acf_resid(m2)[2]

  # Full model with k = 9
  m3 <- bam(
    z_y ~
      s(closest_angle, k = 9) +
      s(closest_angle, by = vowel, k = 9) +
      s(participant, bs = "re") +
      s(file, bs = "re") +
      s(closest_angle, participant, bs = "fs", k = 9),
    AR.start = data$AR.start,
    rho = rho2,
    data = data,
    method = "fREML",
    discrete = TRUE
  )
  rho3 <- acf_resid(m3)[2]

  # Final model
  final <- bam(
    z_y ~
      s(closest_angle, k = 9) +
      s(closest_angle, by = vowel, k = 9) +
      s(participant, bs = "re") +
      s(file, bs = "re") +
      s(closest_angle, participant, bs = "fs", k = 9),
    AR.start = data$AR.start,
    rho = rho3,
    data = data,
    method = "fREML",
    discrete = TRUE
  )

  list(model = final, data = data)
}

sp_ie <- prepare_pair(sp, c("i", "e"), "e")
sp_uo <- prepare_pair(sp, c("u", "o"), "o")

fit_ie <- fit_pair(sp_ie)
fit_uo <- fit_pair(sp_uo)

gamm_spie_final <- fit_ie$model
gamm_spuo_final <- fit_uo$model

summary(gamm_spie_final)
k.check(gamm_spie_final)

summary(gamm_spuo_final)
k.check(gamm_spuo_final)

# -------------------------------------------------------------------------
# 5. GAMM figure
# -------------------------------------------------------------------------

draw_gamm_plots <- function() {
  layout(
    matrix(c(1, 3, 2, 4), nrow = 2, byrow = TRUE),
    heights = c(3, 1)
  )

  par(mar = c(4, 4, 3, 1))
  plot_smooth(
    gamm_spie_final,
    view = "closest_angle",
    plot_all = "vowel",
    cond = list(vowel = c("e", "i")),
    rm.ranef = TRUE,
    shade = TRUE,
    rug = TRUE,
    col = c("e" = "#FA5252", "i" = "#A869C7"),
    main = "/i/ vs /e/",
    xlab = "Angle",
    ylab = "Tongue height",
    sim.ci = TRUE,
    sim.ci.type = "simultaneous"
  )

  par(mar = c(4, 4, 1, 1))
  plot_diff(
    gamm_spie_final,
    view = "closest_angle",
    comp = list(vowel = c("e", "i")),
    rm.ranef = TRUE,
    shade = TRUE,
    shade.col = "grey70",
    rug = TRUE,
    col = "black",
    col.diff = "#D90429",
    xlab = "Angle",
    ylab = "Difference",
    sim.ci = TRUE,
    sim.ci.type = "simultaneous"
  )

  par(mar = c(4, 4, 3, 1))
  plot_smooth(
    gamm_spuo_final,
    view = "closest_angle",
    plot_all = "vowel",
    cond = list(vowel = c("o", "u")),
    rm.ranef = TRUE,
    shade = TRUE,
    rug = TRUE,
    col = c("o" = "#4A9663", "u" = "#2B67E0"),
    main = "/u/ vs /o/",
    xlab = "Angle",
    ylab = "Tongue height",
    sim.ci = TRUE,
    sim.ci.type = "simultaneous"
  )

  par(mar = c(4, 4, 1, 1))
  plot_diff(
    gamm_spuo_final,
    view = "closest_angle",
    comp = list(vowel = c("o", "u")),
    rm.ranef = TRUE,
    shade = TRUE,
    shade.col = "grey70",
    rug = TRUE,
    col = "black",
    col.diff = "#D90429",
    xlab = "Angle",
    ylab = "Difference",
    sim.ci = TRUE,
    sim.ci.type = "simultaneous"
  )
}

png(
  "spanish_gamm_results.png",
  width = 12, height = 6, units = "in", res = 300, bg = "white"
)
draw_gamm_plots()
dev.off()

# Display plot
draw_gamm_plots()