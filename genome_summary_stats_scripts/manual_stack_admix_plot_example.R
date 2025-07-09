### chatGPT's advice for always stacking the largest proportion first.
# geom_rect seems to work

library(ggplot2)
library(dplyr)

# Example data
data <- tibble(
  pop = c("Mallee", "Mallee", "Mallee", "Mallee", "Mallee", "Mallee", "Mallee"),
  Individual = c("Ind16", "Ind16", "Ind16", "Ind14", "Ind14", "Ind14", "Ind18"),
  group_id = c("1", "1", "1", "1", "1", "1", "1"),
  x_pos = c(1, 1, 1, 2, 2, 2, 3),
  population = c("q1", "q2", "q3", "q1", "q2", "q3", "q1"),
  proportion = c(0.3506627, 0.5955527, 0.0537847, 0.3503368, 0.5851296, 0.0645335, 0.3484213)
)

# Step 1: Ensure data is sorted in the correct stacking order
data <- data %>%
  group_by(x_pos) %>%
  arrange(desc(proportion), .by_group = TRUE) %>%
  ungroup()

# Step 2: Create a new column for localized factor levels
data <- data %>%
  group_by(x_pos) %>%
  mutate(population = factor(population, levels = population)) %>%  # Local factor levels per x_pos
  ungroup()

# Step 3: Plot the stacked bar chart
ggplot(data, aes(x = factor(x_pos), y = proportion, fill = population)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Stacked Barplot with Correct Stacking Order",
    x = "Individual Position",
    y = "Proportion",
    fill = "Population"
  ) +
  theme_minimal()


####################################################################################################


library(ggplot2)
library(dplyr)

# Example data
data <- tibble(
  pop = c("Mallee", "Mallee", "Mallee", "Mallee", "Mallee", "Mallee", "Mallee"),
  Individual = c("Ind16", "Ind16", "Ind16", "Ind14", "Ind14", "Ind14", "Ind18"),
  group_id = c("1", "1", "1", "1", "1", "1", "1"),
  x_pos = c(1, 1, 1, 2, 2, 2, 3),
  population = c("q1", "q2", "q3", "q1", "q2", "q3", "q1"),
  proportion = c(0.3506627, 0.5955527, 0.0537847, 0.3503368, 0.5851296, 0.0645335, 0.3484213)
)

# Step 1: Arrange the data for correct stacking order
data <- data %>%
  group_by(x_pos) %>%
  arrange(desc(proportion), .by_group = TRUE) %>%
  mutate(
    ymin = c(0, cumsum(proportion)[-n()]),  # Bottom of the segment
    ymax = cumsum(proportion)              # Top of the segment
  ) %>%
  ungroup()

# Step 2: Plot manually using geom_rect()
ggplot(data) +
  geom_rect(aes(
    xmin = as.numeric(x_pos) - 0.4,   # Bar width left
    xmax = as.numeric(x_pos) + 0.4,   # Bar width right
    ymin = ymin,                      # Bottom of segment
    ymax = ymax,                      # Top of segment
    fill = population
  )) +
  scale_x_continuous(breaks = unique(data$x_pos), labels = unique(data$x_pos)) +
  labs(
    title = "Stacked Barplot with Correct Stacking Order",
    x = "Individual Position",
    y = "Proportion",
    fill = "Population"
  ) +
  theme_minimal()