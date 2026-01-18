####----加载包----####
library(tidyverse)
library(readxl)
library(ggplot2)
library(maps)  # 提供 map_data("world")

####----加载数据----####
# 1️⃣ 读取世界地图数据
world <- map_data("world") %>% 
  filter(region != "Antarctica")

# 2️⃣ 读取你的样本信息
Sample_Info <- read_csv("/Users/ss/Desktop/clean_self_ref_paper_excluded_age.csv", col_names = TRUE)

# 3️⃣ 汇总每个国家的被试人数
country_summary <- Sample_Info %>% 
  group_by(Data_Collection_country) %>% 
  summarise(Total_Sample = sum(Sample_patients, na.rm = TRUE))

# 4️⃣ 合并地图与汇总后的数据
world_with_data <- world %>% 
  left_join(country_summary, by = c("region" = "Data_Collection_country"))

####----绘图----####
p <- ggplot() +
  # 绘制国家边界并填充样本量
  geom_polygon(
    data = world_with_data,
    aes(x = long, y = lat, group = group, fill = Total_Sample),
    color = "gray50", linewidth = 0.3
  ) +
  # 样本量颜色渐变
  scale_fill_gradient(
    low = "#cfe7c4",  # 浅色（样本量少）
    high = "#4f845c", # 深色（样本量多）
    na.value = "white",
    name = "Sample Size"
  ) +
  # 投影方式
  coord_quickmap() +
  # 主题和风格
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.position = "bottom",          # 放到底部
    legend.justification = "center",     # 居中
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  )

####----保存图片----####
ggsave(
  filename = "/Users/ss/Desktop/patients_sample_map.pdf",
  plot = p,
  height = 6,
  width = 12
)

ggsave(
  filename = "/Users/ss/Desktop/patients_sample_map.svg",
  plot = p,
  height = 6,
  width = 12
)

# 展示结果
print(p)

# ================================================================================================
# 计算"disease_name_G2"列中各值的出现次数
counts <- table(df$Disease_name_G2)

# 查看结果
print(counts)

counts <- table(df$First_Author)

# 查看结果
print(counts)

total <- sum(counts)
print(total)

####---- 加载必要包 ----####
library(tidyverse)
library(scales)

####---- 读取数据 ----####
df <- read_csv("/Users/ss/Desktop/clean_self_ref_paper_excluded_age.csv")
df <- df %>% filter(!is.na(First_Author))

####---- 按疾病出现频率排序 ----####
df <- df %>% mutate(Disease_name_G2 = fct_infreq(Disease_name_G2))

####---- 绘图：Bar+散点图（单纵轴，y轴延伸至60） ----####
p <- ggplot(df, aes(x = Disease_name_G2)) +
  geom_bar(
    aes(y = ..count..),
    fill = "gray70",
    color = "black",
    linewidth = 0.3,
    width = 0.5,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(y = Sample_patients, color = Disease_Diagnostic_Criteria_G2),
    position = position_jitterdodge(
      dodge.width = 0.5,
      jitter.width = 0.15,
      jitter.height = 0
    ),
    size = 2.5,
    alpha = 0.7
  ) +
  labs(
    x = "Disorder type",
    y = "Number of psychiatric patients",
    color = "Diagnostic criteria"
  ) +
  # ✅ 在这里新增 sec.axis，右边轴线才会出现
  scale_y_continuous(
    limits = c(0, 60),
    expand = c(0, 0),
    breaks = breaks_extended(n = 6),
    sec.axis = dup_axis(name = "Number of studies")  # <-- 新增
  ) +
  scale_color_brewer(palette = "Paired") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    
    # ✅ 保证右轴线和刻度样式生效
    axis.line.y.right = element_line(color = "black"),
    axis.ticks.y.right = element_line(color = "black"),
    axis.text.y.right = element_text(size = 10, color = "black"),
    axis.title.y.right = element_text(size = 12, face = "bold", color = "black", margin = margin(l = 8)),
    
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.background = element_rect(color = "black", fill = NA, linewidth = 0.2),
    legend.key = element_blank()
  )

print(p)


####----保存图片----####
ggsave(
  filename = "/Users/ss/Desktop/bar_plot.pdf",
  plot = p,
  height = 6,
  width = 12
)

ggsave(
  filename = "/Users/ss/Desktop/bar_plot.svg",
  plot = p,
  height = 6,
  width = 12
)


