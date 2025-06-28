library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
install.packages("svglite")
# 数据准备 ----------------------------------------------------------------
df <- read.csv("C:/Users/自动挡赛车手/Desktop/doublesex文章/vc.csv", 
               sep = ",", 
               header = TRUE,
               check.names = FALSE,
               na.strings = c("NA", ""))

# 预处理步骤 --------------------------------------------------------------
highlight_genes <- na.omit(df$hightlight)

df_clean <- df %>%
  filter(
    !is.na(`diffexp_log2fc_WT_Female-vs-Csdsx_ko_F`),
    !is.na(`diffexp_deseq2_pvalue_WT_Female-vs-Csdsx_ko_F`)
  ) %>%
  mutate(
    inverted_log2fc = `diffexp_log2fc_WT_Female-vs-Csdsx_ko_F`,
    logpvalue = -log10(`diffexp_deseq2_pvalue_WT_Female-vs-Csdsx_ko_F`),
    group = case_when(
      inverted_log2fc > 1 & `diffexp_deseq2_pvalue_WT_Female-vs-Csdsx_ko_F` < 0.05 ~ "up",
      inverted_log2fc < -1 & `diffexp_deseq2_pvalue_WT_Female-vs-Csdsx_ko_F` < 0.05 ~ "down",
      TRUE ~ "none"
    ),
    label = ifelse(
      gene_id %in% highlight_genes,
      coalesce(desc_pfam, gene_id),
      ""
    )
  )

# 颜色映射 --------------------------------------------------------------
df_color <- df_clean %>%
  mutate(
    color = case_when(
      group == "none" & label == "" ~ "color1",
      group == "up"   & label == "" ~ "color3",   # 上调基因基础颜色
      group == "down" & label == "" ~ "color2",   # 下调基因基础颜色
      group == "up"   & label != "" ~ "color5",   # 高亮上调基因
      group == "down" & label != "" ~ "color4",   # 高亮下调基因
      TRUE ~ "color1"
    )
  ) %>%
  arrange(color)

# 可视化绘制 ------------------------------------------------------------
base_plot <- ggscatter(
  df_color,
  x = "inverted_log2fc",
  y = "logpvalue",
  color = "color",
  palette = c("#bcbcbc", "#8abddc", "#ffab84", "#0051a6", "#be000e"), # 调整颜色顺序
  size = 2,
  alpha = 0.7
) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "grey50",
             linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", 
             color = "grey50",
             linewidth = 0.3) +
  labs(
    x = expression(paste(Log[2], " Fold Change ")),
    y = expression(paste(-Log[10], " P-value"))
  ) +
  theme_bw() +
  theme(
    plot.title    = element_blank(),
    panel.grid    = element_blank(),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5), # 完整四周框线
    legend.position = "none", # 移除颜色图例
    axis.ticks    = element_line(color = "black", linewidth = 0.3)
  )

# 智能标签布局（左标签左对齐，右标签右对齐）---------------------------------
final_plot <- base_plot +
  geom_text_repel(
    data = subset(df_color, label != "" & inverted_log2fc < 0),
    aes(label = label),
    color = "black",
    size = 3,
    direction = "y",  
    nudge_x = -20 - subset(df_color, label != "" & inverted_log2fc < 0)$inverted_log2fc, # 左侧标签左推
    segment.color = "black",
    segment.size  = 0.25,
    box.padding   = 0.3,
    max.time      = 10
  ) +
  geom_text_repel(
    data = subset(df_color, label != "" & inverted_log2fc >= 0),
    aes(label = label),
    color = "black",
    size = 3,
    direction = "y",
    nudge_x = 15 - subset(df_color, label != "" & inverted_log2fc >= 0)$inverted_log2fc, # 右侧标签右推
    segment.color = "black",
    segment.size  = 0.25,
    box.padding   = 0.3,
    max.time      = 10
  )
num_up   <- sum(df_clean$group == "up",   na.rm = TRUE)
num_down <- sum(df_clean$group == "down", na.rm = TRUE)
message("显著上调基因数量: ",   num_up)
message("显著下调基因数量: ",   num_down)
message("总差异基因数量: ",     num_up + num_down)
# 保存输出（宽度10cm, 高度12cm）------------------------------------------------
# 保存为 TIFF
ggsave("final_volcano_pro.tiff", 
       plot   = final_plot,
       device = "tiff",
       dpi    = 1200,
       width  = 10,   # 宽度 10 cm
       height = 12,
       units  = "cm")

# 保存为 SVG
ggsave("final_volcano_pro.svg", 
       plot   = final_plot,
       device = "svg",
       width  = 10,   # 宽度 10 cm
       height = 12,
       units  = "cm")