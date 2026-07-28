# 强制跳过所有交互提示
options(warn = -1)
options(menu.graphics = FALSE)

# 清除所有变量
rm(list = ls())

# 静默安装包
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# 安装并加载包
cat("正在安装和加载必要的包...\n")
required_packages <- c("ggplot2", "ggspatial", "sf", "rnaturalearth", "maps", "mapdata")
for (pkg in required_packages) {
  suppressWarnings(suppressMessages(install_if_missing(pkg)))
}
cat("包加载完成！\n\n")

# 创建数据
data <- data.frame(
  SITE = c("SY365BB", "SY366YB", "SY366YW", "SY368YW", "SY456YB", "SY457YW", 
           "SY459YW", "S1_S2S4S3", "SQ_58", "SQ_581", "C1", "C2", "C3", 
           "NS", "S", "ES", "R2111N_S500", "R2111N_S300"),
  Lon = c(111.11707, 111.12198, 111.12114, 111.114929, 111.12076, 111.12076,
          111.05574, 110.4, 110.47, 110.4, 110.4724, 110.4719, 110.4718,
          110.43, 110.5, 110.4, 111, 111),
  Lat = c(17.718045, 17.7036, 17.70335, 17.699036, 17.70335, 17.70335,
          17.62273, 16.9, 16.73, 16.69, 16.7283, 16.7285, 16.7284,
          16.71, 16.5, 16.9, 17.2, 17)
)

# 显示数据信息
cat("=== 数据信息 ===\n")
cat("站点数量:", nrow(data), "\n")
cat("经度范围:", sprintf("%.4f ~ %.4f", min(data$Lon), max(data$Lon)), "\n")
cat("纬度范围:", sprintf("%.4f ~ %.4f", min(data$Lat), max(data$Lat)), "\n\n")

# 设置地图范围（稍微扩大）
lon_range <- c(110.0, 111.5)
lat_range <- c(16.3, 17.8)

# 方法1：使用maps包创建简单地图（避免坐标系统问题）
cat("正在生成地图...\n")

# 创建输出目录
output_dir <- "~/Desktop/R_maps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 创建基础地图（使用maps包，更简单）
# 获取中国地图数据
china_map <- map_data("world", region = "China")

# 创建地图 - 简单版本
map_simple <- ggplot() +
  # 绘制中国地图
  geom_polygon(data = china_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "gray60", linewidth = 0.3) +
  
  # 绘制海洋背景
  geom_rect(aes(xmin = lon_range[1], xmax = lon_range[2],
                ymin = lat_range[1], ymax = lat_range[2]),
            fill = "lightblue", alpha = 0.2) +
  
  # 绘制站点
  geom_point(data = data, aes(x = Lon, y = Lat),
             color = "red", size = 3, alpha = 0.8) +
  
  # 添加站点标签
  geom_text(data = data, aes(x = Lon, y = Lat, label = SITE),
            hjust = 0, vjust = 1.5, size = 2.5, fontface = "bold",
            check_overlap = FALSE) +
  
  # 设置坐标范围
  coord_fixed(ratio = 1.5, 
              xlim = lon_range, 
              ylim = lat_range,
              expand = TRUE) +
  
  # 添加比例尺和指北针
  annotation_scale(location = "br", width_hint = 0.3) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering) +
  
  # 标题和标签
  labs(title = "站点位置分布图",
       subtitle = paste("共", nrow(data), "个站点"),
       x = "经度 (°E)",
       y = "纬度 (°N)") +
  
  # 主题
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.background = element_rect(fill = "lightblue", alpha = 0.1),
    panel.grid = element_line(color = "gray90", linewidth = 0.3)
  )

# 显示地图
print(map_simple)

# 保存为PDF
pdf_path1 <- file.path(output_dir, "sites_map_simple.pdf")
ggsave(pdf_path1, plot = map_simple, width = 10, height = 8, device = "pdf")
cat("✓ 简单地图已保存:", pdf_path1, "\n")

# 方法2：使用ggplot2 + 海岸线数据
cat("正在生成海岸线地图...\n")

# 获取世界海岸线数据
world_coast <- map_data("world")

# 创建海岸线地图
map_coast <- ggplot() +
  # 绘制海洋背景
  geom_rect(aes(xmin = lon_range[1], xmax = lon_range[2],
                ymin = lat_range[1], ymax = lat_range[2]),
            fill = "lightblue", alpha = 0.3) +
  
  # 绘制海岸线
  geom_polygon(data = world_coast,
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "gray50", linewidth = 0.3) +
  
  # 绘制站点（按经度大小着色）
  geom_point(data = data, aes(x = Lon, y = Lat, color = Lon),
             size = 4, alpha = 0.8) +
  
  # 添加站点标签
  geom_text(data = data, aes(x = Lon, y = Lat, label = SITE),
            hjust = 0, vjust = 1.5, size = 2.5, fontface = "bold") +
  
  # 颜色渐变
  scale_color_gradient(low = "blue", high = "red", name = "经度") +
  
  # 设置坐标范围
  coord_fixed(ratio = 1.5, 
              xlim = lon_range, 
              ylim = lat_range,
              expand = TRUE) +
  
  # 添加比例尺和指北针
  annotation_scale(location = "br", width_hint = 0.3) +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering) +
  
  # 标题
  labs(title = "站点位置分布图（海岸线）",
       subtitle = paste("共", nrow(data), "个站点"),
       x = "经度 (°E)",
       y = "纬度 (°N)") +
  
  # 主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "lightblue", alpha = 0.1)
  )

# 保存海岸线地图
pdf_path2 <- file.path(output_dir, "sites_map_coast.pdf")
ggsave(pdf_path2, plot = map_coast, width = 10, height = 8, device = "pdf")
cat("✓ 海岸线地图已保存:", pdf_path2, "\n")

# 方法3：更详细的区域地图
cat("正在生成详细区域地图...\n")

# 获取更详细的中国省份数据
if (requireNamespace("maps", quietly = TRUE)) {
  china_province <- map_data("china")
  
  map_detail <- ggplot() +
    # 海洋背景
    geom_rect(aes(xmin = lon_range[1], xmax = lon_range[2],
                  ymin = lat_range[1], ymax = lat_range[2]),
              fill = "#E6F3FF", alpha = 0.5) +
    
    # 中国省份
    geom_polygon(data = china_province,
                 aes(x = long, y = lat, group = group),
                 fill = "#F5F5DC", color = "gray40", linewidth = 0.3) +
    
    # 站点（按分组着色）
    # 创建分组
    data$group <- ifelse(data$Lon > 110.8, "东部", "西部")
  
  geom_point(data = data, aes(x = Lon, y = Lat, fill = group),
             size = 4, shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
    
    # 站点标签
    geom_text(data = data, aes(x = Lon, y = Lat, label = SITE),
              hjust = 0, vjust = 1.5, size = 2.5, fontface = "bold") +
    
    # 填充颜色
    scale_fill_manual(values = c("东部" = "red", "西部" = "blue"), name = "区域") +
    
    # 设置坐标范围
    coord_fixed(ratio = 1.5, 
                xlim = lon_range, 
                ylim = lat_range,
                expand = TRUE) +
    
    # 添加比例尺和指北针
    annotation_scale(location = "br", width_hint = 0.3) +
    annotation_north_arrow(location = "tr",
                           style = north_arrow_fancy_orienteering) +
    
    # 标题
    labs(title = "站点位置分布图（详细区域）",
         subtitle = paste("共", nrow(data), "个站点"),
         x = "经度 (°E)",
         y = "纬度 (°N)") +
    
    # 主题
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.background = element_rect(fill = "#E6F3FF"),
      legend.position = "bottom"
    )
  
  # 保存详细地图
  pdf_path3 <- file.path(output_dir, "sites_map_detail.pdf")
  ggsave(pdf_path3, plot = map_detail, width = 10, height = 8, device = "pdf")
  cat("✓ 详细区域地图已保存:", pdf_path3, "\n")
  
  # 显示详细地图
  print(map_detail)
}

# 保存数据为CSV
csv_path <- file.path(output_dir, "sites_data.csv")
write.csv(data, csv_path, row.names = FALSE)
cat("✓ 数据文件已保存:", csv_path, "\n")

# 创建数据摘要文件
summary_path <- file.path(output_dir, "data_summary.txt")
sink(summary_path)
cat("=== 站点数据统计 ===\n\n")
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("站点总数:", nrow(data), "\n\n")

cat("经纬度范围:\n")
cat("  经度: ", min(data$Lon), "~", max(data$Lon), "\n")
cat("  纬度: ", min(data$Lat), "~", max(data$Lat), "\n\n")

cat("站点列表:\n")
for (i in 1:nrow(data)) {
  cat(sprintf("  %2d. %-12s 经度: %10.5f  纬度: %10.5f\n",
              i, data$SITE[i], data$Lon[i], data$Lat[i]))
}

cat("\n区域分布:\n")
if (exists("data$group")) {
  cat("  东部区域:", sum(data$group == "东部"), "个站点\n")
  cat("  西部区域:", sum(data$group == "西部"), "个站点\n")
}
sink()

cat("✓ 数据摘要已保存:", summary_path, "\n")

# 显示最终汇总信息
cat("\n", strrep("=", 60), "\n")
cat("地图生成完成！\n")
cat(strrep("=", 60), "\n\n")

cat("生成的文件清单:\n")
cat("1. 简单地图:        ", pdf_path1, "\n")
cat("2. 海岸线地图:      ", pdf_path2, "\n")
if (file.exists(pdf_path3)) {
  cat("3. 详细区域地图:    ", pdf_path3, "\n")
}
cat("4. 数据文件 (CSV):  ", csv_path, "\n")
cat("5. 数据摘要 (TXT):  ", summary_path, "\n\n")

cat("数据统计:\n")
cat("- 站点总数:      ", nrow(data), "\n")
cat("- 经度平均值:    ", round(mean(data$Lon), 4), "\n")
cat("- 纬度平均值:    ", round(mean(data$Lat), 4), "\n")
cat("- 中心位置:      (", round(mean(data$Lon), 4), ", ", 
    round(mean(data$Lat), 4), ")\n\n", sep = "")

cat("✓ 所有地图和文件已成功生成！\n")
cat("✓ 请查看桌面上的 R_maps 文件夹\n")
cat("✓ 地图已显示在屏幕上\n")