# 安装并加载必要的包
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("forcats", quietly = TRUE)) install.packages("forcats")
library(ggplot2)
library(dplyr)
library(forcats)

# 1. 设置工作目录
setwd("C:/Users/自动挡赛车手/Desktop/doublesex文章/染色体log2fc")

# 2. 读入三张表
chr_len  <- read.table("cs.txt", sep = "\t", header = FALSE,
                       col.names = c("Chr", "Length"),
                       stringsAsFactors = FALSE)
gene_pos <- read.table("gene_pos.txt", sep = "\t", header = FALSE,
                       col.names = c("Gene", "Chr", "Start", "End"),
                       stringsAsFactors = FALSE)
gene_col <- read.table("gene_col1.txt", sep = "\t", header = FALSE,
                       col.names = c("Gene", "RGB"),
                       stringsAsFactors = FALSE)

# 3. 去除 “ou” 前缀，让 Chr 只剩 “chr…” 
#    否则合并后会出现匹配不到的问题
chr_len$Chr   <- sub("^ou", "", chr_len$Chr, ignore.case = TRUE)
gene_pos$Chr  <- sub("^ou", "", gene_pos$Chr, ignore.case = TRUE)

# 4. 转换颜色：R,G,B -> hex
gene_col$RGB  <- gsub('[" ]', '', gene_col$RGB)
rgb_list      <- strsplit(gene_col$RGB, ",")
gene_col$hex  <- sapply(rgb_list, function(v) {
  grDevices::rgb(as.numeric(v[1]),
                 as.numeric(v[2]),
                 as.numeric(v[3]),
                 maxColorValue = 255)
})

# 5. 合并位置与颜色
dat <- merge(gene_pos, gene_col[, c("Gene", "hex")], by = "Gene")

# 6. 根据 hex 映射到三组，并更新分组标签
color_map <- c(
  "#00ABB8" = "Down-regulated (<-2)",
  "#FFFEFE" = "None (-2 to 2)",
  "#FF9027" = "Up-regulated (>2)"
)
dat$Group <- factor(color_map[toupper(dat$hex)],
                    levels = c("Down-regulated (<-2)",
                               "None (-2 to 2)",
                               "Up-regulated (>2)"))

# 7. 按染色体数字排序并生成新标签
get_chr_num <- function(chr_names) {
  as.numeric(sub("^chr", "", chr_names, ignore.case = TRUE))
}
chr_len <- chr_len %>%
  mutate(Chr_Num = get_chr_num(Chr)) %>%
  arrange(Chr_Num) %>%
  select(-Chr_Num)

# 7.5 为染色体生成新的标签（chr1–chr30），并将 chr895/chr923 重命名为 chrZ/chrW
chr_len <- chr_len %>%
  mutate(Label = case_when(
    Chr == "chr895" ~ "chrZ",
    Chr == "chr923" ~ "chrW",
    TRUE ~ paste0("chr", row_number())
  ))

# 8. 按两行15列布局分配位置
chr_len <- chr_len %>%
  mutate(Row = ceiling(row_number() / 15),    # 每行 15 条
         Col = (row_number() - 1) %% 15 + 1)

# 9. 为基因数据添加分组信息（与染色体对应）
dat <- dat %>%
  left_join(chr_len %>% select(Chr, Row, Col), by = "Chr")

# 10. 计算全局最大长度与标签位置
global_max_length <- max(chr_len$Length)
chr_labels <- chr_len %>%
  group_by(Row) %>%
  mutate(Max_Length = global_max_length) %>%
  ungroup()

# 11. 绘图：柱体更窄、底部留白更大，让“chr”标签完整显示
p <- ggplot() +
  # 染色体柱子（宽度：从 0.15 → 0.10）
  geom_rect(data = chr_len,
            aes(xmin = Col - 0.10, xmax = Col + 0.10,
                ymin = 0, ymax = Length),
            fill = "gray90", color = "gray30", size = 0.2) +
  
  # 基因色块（宽度：从 0.14 → 0.10）
  geom_rect(data = dat,
            aes(xmin = Col - 0.10, xmax = Col + 0.10,
                ymin = Start, ymax = End,
                fill = Group),
            color = NA, alpha = 0.9) +
  
  # 染色体底部标签：新 Label（chr1, chr2…, chrZ, chrW）
  geom_text(data = chr_labels,
            aes(x = Col, y = -global_max_length * 0.07, label = Label),
            size = 2, vjust = 1, fontface = "bold") +
  
  # 染色体顶部标签：原本的 Chr 名称
  geom_text(data = chr_labels,
            aes(x = Col, y = global_max_length * 1.04, label = Chr),
            size = 2, vjust = 0, fontface = "plain") +
  
  # 分面：两行对齐，scales="fixed" 保持纵坐标一致
  facet_wrap(~ Row, nrow = 2, scales = "fixed") +
  
  # 颜色与图例
  scale_fill_manual(name = "Expression",
                    values = c(
                      "Down-regulated (<-2)" = "#00ABB8",
                      "None (-2 to 2)"      = "#A0A0A0",
                      "Up-regulated (>2)"   = "#FF9027"
                    )) +
  
  # 保留左侧 y 轴坐标（以 Mb 为单位），并扩大上下留白
  scale_y_continuous(
    name   = "Genomic position (Mb)",
    labels = function(x) sprintf("%.0f", x/1e6),
    breaks = pretty(c(0, global_max_length), n = 5),
    limits = c(-global_max_length * 0.08, global_max_length * 1.06),
    expand = c(0, 0)
  ) +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL) +
  
  # 主题：紧凑版，同时给底部留出更多空间
  theme_minimal(base_size = 8) +
  theme(
    panel.grid       = element_blank(),
    panel.spacing    = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "top",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.5, "cm"),
    # 加大下边距：保证 “chr” 标签不会被裁剪
    plot.margin      = margin(1.0, 0.5, 1.0, 1.0, "cm")
  )

# 12. 预览与保存
ggsave("chromosome_plot_shrunk_labels_ok.tiff", p,
       width = 6, height = 5, dpi = 600, device = "tiff")
print(p)
