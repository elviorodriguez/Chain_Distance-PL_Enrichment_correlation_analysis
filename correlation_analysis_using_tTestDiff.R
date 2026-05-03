library(tidyverse)

################################################################################
######################### Compute euclidean differences ########################
################################################################################

condition_to_id_correspondance <- c(
  "C4B63_57g96" = "BDF2",
  "C4B63_40g26" = "BDF3",
  "C4B63_25g318" = "BDF4",
  "C4B63_1g51" = "BDF5",
  "BCY84_19389" = "BDF6",
  "C4B63_7g80" = "BDF8"
)

# Load data
enrichments <- read_delim("../bdfs_enrichment/interactors_and_proximal/network.txt", delim = "\t") %>% 
  select(SOURCE, TARGET, `N: Student's T-test Difference Bait_GFP`, `N: Student's T-test Difference Bait_NLS`) %>% 
  rename(diff_GFP = `N: Student's T-test Difference Bait_GFP`,
         diff_NLS = `N: Student's T-test Difference Bait_NLS`) %>% 
  mutate(euclidean_diff = sqrt(diff_GFP^2 + diff_NLS^2),
         BAIT = condition_to_id_correspondance[SOURCE])

################################################################################
######################## Read protein distances tables #########################
################################################################################

# Read distance data (AF3)
distances_table <-  read_delim("./results_af3/distances_averages_20251124_132357.tsv", delim = "\t")

# Reorder SOURCE and TARGET
distances_table_reordered <- distances_table %>%
  mutate(
    source_in_conditions = SOURCE %in% names(condition_to_id_correspondance),
    target_in_conditions = TARGET %in% names(condition_to_id_correspondance)
  ) %>%
  # Case 1: Both are in conditions - duplicate the row
  bind_rows(
    filter(., source_in_conditions & target_in_conditions)
  ) %>%
  # Switch SOURCE and TARGET when needed
  mutate(
    needs_switch = !source_in_conditions & target_in_conditions,
    SOURCE_new = if_else(needs_switch, TARGET, SOURCE),
    TARGET_new = if_else(needs_switch, SOURCE, TARGET),
    seq1_new = if_else(needs_switch, seq2, seq1),
    seq2_new = if_else(needs_switch, seq1, seq2)
  ) %>%
  # For duplicated rows (both in conditions), switch one of them
  group_by(SOURCE, TARGET) %>%
  mutate(
    row_num = row_number(),
    should_switch_duplicate = source_in_conditions & target_in_conditions & row_num == 2,
    SOURCE_new = if_else(should_switch_duplicate, TARGET, SOURCE_new),
    TARGET_new = if_else(should_switch_duplicate, SOURCE, TARGET_new),
    seq1_new = if_else(should_switch_duplicate, seq2, seq1_new),
    seq2_new = if_else(should_switch_duplicate, seq1, seq2_new)
  ) %>%
  ungroup() %>%
  # Keep only the new columns
  select(SOURCE = SOURCE_new, TARGET = TARGET_new, seq1 = seq1_new, seq2 = seq2_new,
         mean_distance, SD_distance, median_distance, IQR_distance)

# View the result
distances_table_reordered

# Merge the datasets by SOURCE and TARGET
merged_data <- enrichments %>%
  inner_join(distances_table_reordered, by = c("SOURCE", "TARGET"))

# Check the merged data
print(merged_data)
print(paste("Number of matched pairs:", nrow(merged_data)))

################################################################################
######################## Compute and plot correlation ##########################
################################################################################

# Function to create correlation plot
create_correlation_plot <- function(data, x_var, y_var, color_var, x_label, y_label,  color_label) {
  cor_result <- cor.test(data[[x_var]], data[[y_var]], method = "spearman", exact = FALSE)
  rho <- round(cor_result$estimate, 3)
  p_value <- round(cor_result$p.value, 3)
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[color_var]])) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +
    geom_point(aes(fill = .data[[color_var]]), shape = 21, color = "black", stroke = 0.8, alpha = 0.8, size = 4) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("ρ = ", rho, " \np-val = ", ifelse(p_value < 0.001, "< 0.001", p_value), " "),
             hjust = 1, vjust = 1.1, size = 6, color = "blue", fontface = "bold") +
    theme_bw() +
    labs(x = x_label, y = y_label, color = color_label) +
    theme(
      text = element_text(size = 24),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )
  
}

# Create the correlation plot
correlation_plot <- create_correlation_plot(
  data = merged_data,
  x_var = "mean_distance",
  y_var = "diff_GFP",
  color_var = "BAIT",
  x_label = "Mean Center of Mass Distance (Å)",
  y_label = "Enrichment (T-test Difference)",
  color_label = "Bait"
)

# Display the plot
print(correlation_plot)

# Optional: Save the plot
ggsave("correlation_distance_vs_enrichment.png",
       correlation_plot,
       width = 7,
       height = 6,
       dpi = 300)

# Optional: Print detailed correlation statistics
cor_result <- cor.test(merged_data$mean_distance, 
                       merged_data$diff_GFP, 
                       method = "spearman", 
                       exact = FALSE)
print(cor_result)

# Optional: Also check Pearson correlation
cor_result_pearson <- cor.test(merged_data$mean_distance, 
                               merged_data$diff_GFP, 
                               method = "pearson")
print(cor_result_pearson)


################################################################################
############### Predict the distance between core and HAT modules ##############
################################################################################

# Fit model
model <- lm(mean_distance ~ diff_GFP, data = merged_data)
summary(model)

source_names <- c(
  "C4B63_40g26" = "BDF3",
  "C4B63_25g318" = "BDF4",
  "C4B63_1g51" = "BDF5"
  # "C4B63_7g80" = "BDF8"
)

target_names <- c(
  "C4B63_6g222" = "HAT2",
  "C4B63_89g65" = "EPL2ap",
  "C4B63_242g4" = "EPL2"
)

selected_pairs <- enrichments %>% 
  filter(
    SOURCE %in% names(source_names),
    TARGET %in% names(target_names)
  ) %>%
  mutate(
    source_name = source_names[SOURCE],
    target_name = target_names[TARGET]
  )


selected_pairs <- selected_pairs %>%
  mutate(predicted_distance = predict(model, newdata = selected_pairs))


final_table <- selected_pairs %>%
  select(SOURCE, TARGET, source_name, target_name, diff_GFP, predicted_distance)



write.csv(final_table, "predicted_distances.csv", row.names = FALSE)

#install.packages("igraph")
library(igraph)
#install.packages("ggraph")
library(ggraph)

g <- graph_from_data_frame(
  final_table,
  directed = FALSE
)

E(g)$weight <- 1 / final_table$predicted_distance
layout_coords <- layout_with_fr(g, weights = E(g)$weight)
layout_df <- as.data.frame(layout_coords)
colnames(layout_df) <- c("x", "y")
layout_df$protein_id <- V(g)$name
E(g)$label <- paste0(round(final_table$predicted_distance, 1), " Å")
V(g)$color <- c("purple", "purple", "purple", "pink","pink","pink")

ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
  geom_edge_link(aes(label = label),
                 angle_calc = "along",
                 label_dodge = unit(2, "mm"),
                 label_size = 4,
                 edge_width = 0.6) +
  geom_node_point(
    aes(color = V(g)$color),
    size = 22,
    shape = 21,          # required for fill + stroke
    fill = V(g)$color,   # inner fill color
    stroke = 1,          # border thickness
    colour = "black"     # border color
  ) +
  geom_node_text(aes(label = c(source_names, target_names)[name])) +
  theme_void()+ 
  expand_limits(
    x = c(min(layout_df$x) - 0.5, max(layout_df$x) + 0.5),
    y = c(min(layout_df$y) - 0.5, max(layout_df$y) + 0.5)
  )







# --- assume final_table already exists and has these columns ---
# SOURCE, TARGET, source_name, target_name, diff_GFP, predicted_distance

# --- create graph ---
g <- graph_from_data_frame(final_table, directed = FALSE)

# --- assign colors to nodes (adjust as needed) ---
V(g)$color <- c(rep("purple", 3), rep("pink", 3))  # 3 sources + 3 targets

# --- soften closeness to prevent collapse ---
k <- 0.05  # adjust this to control intensity
E(g)$length <- final_table$predicted_distance / (1 + k * final_table$predicted_distance)

# --- compute pairwise shortest-path distances ---
sp_dist <- distances(g, weights = E(g)$length)

# Smooth spreading with adjustable midpoint
midpoint <- median(sp_dist[is.finite(sp_dist)])
steepness <- 20

sp_dist <- midpoint * (sp_dist / midpoint)^steepness

# handle disconnected graph
if (any(is.infinite(sp_dist))) {
  max_finite <- max(sp_dist[is.finite(sp_dist)])
  sp_dist[is.infinite(sp_dist)] <- max_finite * 2
}

# --- classical MDS layout based on softened distances ---
mds_coords <- cmdscale(as.dist(sp_dist), k = 2, eig = TRUE)
layout_df <- as.data.frame(mds_coords$points)
colnames(layout_df) <- c("x", "y")
layout_df$protein_id <- rownames(layout_df)
layout_df <- layout_df[match(V(g)$name, layout_df$protein_id), ]

# --- add a small margin ---
layout_df <- layout_df %>%
  mutate(
    x = scales::rescale(x, to = c(0.1, 0.9)),
    y = scales::rescale(y, to = c(0.1, 0.9))
  )

# --- compute actual Euclidean lengths (optional) ---
edge_ends <- ends(g, E(g))
u_idx <- match(edge_ends[,1], layout_df$protein_id)
v_idx <- match(edge_ends[,2], layout_df$protein_id)
E(g)$actual_length <- sqrt((layout_df$x[u_idx] - layout_df$x[v_idx])^2 +
                             (layout_df$y[u_idx] - layout_df$y[v_idx])^2)

# --- edge labels ---
E(g)$label <- paste0(round(final_table$predicted_distance, 1), " Å")

# --- plot ---
ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
  geom_edge_link(aes(label = label),
                 angle_calc = "along",
                 label_dodge = unit(2, "mm"),
                 label_size = 4,
                 edge_width = 0.6) +
  geom_node_point(
    aes(color = V(g)$color),
    size = 22,
    shape = 21,
    fill = V(g)$color,
    stroke = 1,
    colour = "black"
  ) +
  geom_node_text(aes(label = c(source_names, target_names)[name]),
                 size = 5) +
  theme_void() +
  expand_limits(
    x = c(min(layout_df$x) - 0.05, max(layout_df$x) + 0.05),
    y = c(min(layout_df$y) - 0.05, max(layout_df$y) + 0.05)
  )
