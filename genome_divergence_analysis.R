####################################################
# ENHANCED VISUALIZATION OF GENOME DIVERGENCE MATRIX #
####################################################

# 0. THIẾT LẬP THƯ MỤC ----
setwd("/Users/shinra/Mash_output")

# 1. KIỂM TRA FILE ----
if (!file.exists("divergence_matrix.tsv")) {
  stop("⛔ Không tìm thấy file divergence_matrix.tsv trong:\n", getwd())
}

# 2. CÀI ĐẶT THƯ VIỆN ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggrepel, viridis, scales, ggridges, patchwork, Hmisc)

# 3. ĐỌC VÀ TIỀN XỬ LÝ DỮ LIỆU ----
div_data <- read_tsv("divergence_matrix.tsv", show_col_types = FALSE) %>% 
  filter(Sample1 != Sample2, !is.na(Distance), is.finite(Distance)) %>%
  # Sửa tên cột chuẩn xác
  mutate(Divergence_Group = cut(Distance,
                                breaks = c(0, 0.03, 0.06, 0.09, 0.12),
                                labels = c("0-3%", "3-6%", "6-9%", "9-12%"))) %>%
  rename(Divergence_Group = Divergence_Group)  # Đảm bảo tên cột nhất quán

# Thêm đoạn kiểm tra
cat("\n✅ KIỂM TRA DIVERGENCE GROUP:\n")
print(table(div_data$Divergence_Group))

# 4. TÍNH TOÁN THỐNG KÊ ----
# Thống kê tổng quan
overall_stats <- div_data %>%
  summarise(
    Mean = mean(Distance),
    SD = sd(Distance),
    Min = min(Distance),
    Max = max(Distance),
    Median = median(Distance),
    Q1 = quantile(Distance, 0.25),
    Q3 = quantile(Distance, 0.75)
  )

# Thống kê theo từng cặp mẫu
pairwise_stats <- div_data %>%
  group_by(Sample1, Sample2) %>%
  summarise(Distance = first(Distance), .groups = 'drop')

# 5. VẼ BIỂU ĐỒ PHÂN BỐ DIVERGENCE ----
dist_plot <- ggplot(div_data, aes(x = Distance)) +
  # Chỉ hiển thị đường density nhưng scale theo count
  geom_density(
    aes(y = after_stat(count)),  # Giữ nguyên đơn vị count
    color = "black",
    linewidth = 1.2,
    fill = NA  # Không tô màu vùng dưới đường
  ) +
  
  # Thống kê mean
  geom_vline(xintercept = overall_stats$Mean, linetype = "dashed", color = "red") +
  annotate("text",
           x = overall_stats$Mean + 0.005,
           y = Inf,
           label = paste("Mean =", round(overall_stats$Mean, 4)),
           vjust = 2,
           color = "red") +
  
  # Định dạng trục
  scale_x_continuous(
    labels = percent_format(),
    limits = c(0, 0.12),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 1)  # Hiển thị số nguyên
  ) +
  
  labs(
    title = "Distribution of Pairwise Divergences",
    subtitle = "Frequency distribution of divergence values",
    x = "Divergence",
    y = "Count of Sample Pairs"  # Giữ nguyên tên trục Y
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )

print(dist_plot)
# 6. VẼ BIỂU ĐỒ MA TRẬN HEATMAP ----
# Chọn một subset các mẫu để tránh quá tải (nếu có nhiều mẫu)
selected_samples <- div_data %>%
  group_by(Sample1) %>%
  summarise(mean_div = mean(Distance)) %>%
  arrange(desc(mean_div)) %>%
  slice_head(n = 30) %>%  # Giới hạn 30 mẫu có divergence cao nhất
  pull(Sample1)

heatmap_data <- div_data %>%
  filter(Sample1 %in% selected_samples, Sample2 %in% selected_samples)

heatmap_plot <- ggplot(heatmap_data, aes(x = Sample1, y = Sample2, fill = Distance)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(0, 0.12), 
                     labels = percent_format()) +
  labs(title = "Pairwise Divergence Heatmap (Top 30 Samples)",
       x = NULL, y = NULL, fill = "Divergence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right")

# 7. VẼ BIỂU ĐỒ XẾP HẠNG DIVERGENCE ----
rank_plot <- div_data %>%
  group_by(Sample1) %>%
  summarise(mean_div = mean(Distance)) %>%
  arrange(desc(mean_div)) %>%
  mutate(Rank = row_number()) %>%
  ggplot(aes(x = Rank, y = mean_div)) +
  geom_point(aes(color = ifelse(mean_div > 0.06, "High", "Low")), size = 2) +
  geom_text_repel(aes(label = ifelse(mean_div > 0.06, Sample1, "")),
                  size = 3, max.overlaps = 20) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = c("High" = "red", "Low" = "blue")) +
  labs(title = "Sample Ranking by Mean Divergence",
       y = "Mean Divergence",
       x = "Rank") +
  theme_minimal() +
  theme(legend.position = "none")

# 8. KẾT HỢP CÁC BIỂU ĐỒ ----
combined_plot <- (dist_plot | rank_plot) / heatmap_plot +
  plot_annotation(title = "Comprehensive Analysis of Genome Divergences",
                  subtitle = paste("Overall statistics: Mean =", round(overall_stats$Mean, 4),
                                   "| SD =", round(overall_stats$SD, 4),
                                   "| Range =", round(overall_stats$Min, 4), "-", round(overall_stats$Max, 4)))

# 9. XUẤT KẾT QUẢ ----
ggsave("divergence_comprehensive_plot.png", combined_plot, 
       width = 16, height = 12, dpi = 300)

# Xuất thống kê
write_tsv(overall_stats, "divergence_overall_stats.tsv")
write_tsv(pairwise_stats, "divergence_pairwise_stats.tsv")

# In thống kê ra console
cat("\n📊 OVERALL STATISTICS:\n")
cat("-----------------------\n")
cat(sprintf("Mean divergence: %.4f\n", overall_stats$Mean))
cat(sprintf("Standard deviation: %.4f\n", overall_stats$SD))
cat(sprintf("Range: %.4f - %.4f\n", overall_stats$Min, overall_stats$Max))
cat(sprintf("Median: %.4f\n", overall_stats$Median))
cat(sprintf("IQR: %.4f - %.4f\n", overall_stats$Q1, overall_stats$Q3))

message("\n✅ HOÀN TẤT!",
        "\n- Biểu đồ tổng hợp: divergence_comprehensive_plot.png",
        "\n- Thống kê tổng quan: divergence_overall_stats.tsv",
        "\n- Thống kê theo cặp: divergence_pairwise_stats.tsv")