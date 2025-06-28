library(pheatmap)

# 1. 正确读取数据（处理数值类型问题）
full_data <- read.delim(
  "C:\\Users\\自动挡赛车手\\Desktop\\doublesex文章\\DSX-KO1.txt", 
  header = TRUE, 
  check.names = FALSE
)

# 设置行名为 Sample 列
rownames(full_data) <- full_data$Sample

# 只保留数值列（假设前两列是 Sample 和 Group）
data <- full_data[, -c(1, 2)]

# 2. 自定义样本顺序
custom_order <- c(
  # dsx-ko-F 样本
  "4.3",  "5.5",  "28.1", "28.3",
  # WT-F 样本
  "WT-F1", "WT-F2", "WT-F3",
  # WT-M 样本
  "WT-M1", "WT-M2", "WT-M3"
)

# 3. 按行顺序重排（行是样本）
data_ordered <- data[custom_order, ]

# 4. 创建行注释（样本分组）
sample_groups <- data.frame(
  Group = factor(
    ifelse(
      rownames(data_ordered) %in% c("4.3", "5.5", "28.1", "28.3"),
      "dsx-ko-F",
      ifelse(
        rownames(data_ordered) %in% c("WT-F1", "WT-F2", "WT-F3"),
        "WT-F",
        "WT-M"
      )
    ),
    levels = c("dsx-ko-F", "WT-F", "WT-M")
  ),
  row.names = rownames(data_ordered)
)

# 5. 以 WT-F 组平均作为参照，计算每列的 fold-change（非 log 版本）
#    首先计算 WT-F 组每个基因的平均值
wtf_mean <- colMeans(data_ordered[c("WT-F1", "WT-F2", "WT-F3"), , drop = FALSE])

#    然后对 data_ordered 的每一列做“除以 wtf_mean”的操作
data_fc <- sweep(data_ordered, 2, wtf_mean, FUN = "/")

# （可选）如果想画 log2(FC)，可以加这一行：
# data_fc <- log2(data_fc)

# 6. 自定义热图配色（红-白-蓝渐变）
custom_colors <- colorRampPalette(
  c(
    rgb(230, 31, 24,   maxColorValue = 255),  # 红
    rgb(253, 245, 247, maxColorValue = 255),  # 白
    rgb(23, 72, 156,   maxColorValue = 255)   # 蓝
  )
)(100)

# 7. 输出到 TIFF 文件
tiff(
  "dsx_ko_heatmap_fc.tiff", 
  width       = 20,     # 单位：cm
  height      = 15, 
  units       = "cm", 
  res         = 1200,   # dpi
  compression = "lzw"
)

# 8. 绘制热图
pheatmap(
  mat                = data_fc,         # 用 fold-change 矩阵
  color              = custom_colors,
  cluster_rows       = FALSE,           # 样本（行）不聚类
  cluster_cols       = FALSE,           # 基因（列）不聚类
  show_rownames      = TRUE,            # 显示行名（样本）
  show_colnames      = TRUE,            # 显示列名（基因）
  annotation_row     = sample_groups,   # 行注释：样本分组
  annotation_colors  = list(
    Group = c(
      "dsx-ko-F" = "gray30",
      "WT-F"     = "#FFB6C1",
      "WT-M"     = "#87CEFA"
    )
  ),
  cellwidth          = 20,              # 单元格宽度
  cellheight         = 12,              # 单元格高度
  gaps_row           = c(4, 7),         # 在第4行（dsx-ko-F后）和第7行（WT-F后）添加分隔线
  fontsize_row       = 8,               # 行名字体大小
  fontsize_col       = 8,               # 列名字体大小
  border_color       = NA,              # 不画边框
  main               = ""               # 不显示标题
)

dev.off()  # 关闭图形设备

