#20250507作图 激素作图#
# data1<-X20250506_RCT_endoforR

library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)
library(coin)
library(readr)
library(readxl)

data1<-read_csv("/data/20250506_RCT_endoforR.csv")
attach(data1)
data1 <- data1 %>%
  mutate(
    time = factor(time),
    group = factor(group)
  )
data1$time <- as.numeric(as.character(data1$time))

#TT
results_tt <- data1 %>%
  group_by(time) %>%
  rstatix::wilcox_test(TT ~ group) %>%  # 明确指定包
  left_join(
    data1 %>% 
      group_by(time) %>% 
      rstatix::wilcox_effsize(TT ~ group),
    by = c("time", "group1", "group2")
  )
print(results_tt)


data1 <- data1 %>% 
  mutate(TT = ifelse(TT > 9, NA, TT))

summary_data <- data1 %>%
  group_by(group, time) %>%
  summarise(
    Mean = mean(TT,na.rm = TRUE),
    SD = sd(TT,na.rm = TRUE),
    .groups = "drop"
  )


tt<-ggplot() +
  geom_line(
    data = summary_data,
    aes(x = time, y = Mean, group = group, color = group),
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    data = summary_data,
    aes(x = time, y = Mean, color = group),
    size = 3,
    position = position_dodge(width = 0.3)
  ) +
  geom_errorbar(
    data = summary_data,
    aes(x = time, ymin = Mean - SD, ymax = Mean + SD, color = group),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  geom_jitter(
    data = data1,
    aes(x = time, y = TT, color = group),
    alpha = 0.4,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.3
    )
  ) +
  labs(
    x = "The diameter of follicles (mm)",
    y = "Serum TT levels (nmol/L)"
  ) +
  scale_color_manual(
    values = c("expri" = "#E41A1C", "control" = "#377EB8"),
    labels = c("expri" = "mLP", "control" = "GnRH-ant")
  ) +
  scale_y_continuous(
    limits = c(-1, 7),
    breaks = seq(0, 7, by = 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(limits = c(0.5,6.5),
                     breaks = c(1,2,3,4,5,6),
                     labels = c("Initial", "<10", "11-12", "13-15", "16-17", "≥18")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),  # 加粗坐标轴线
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),  # 添加 Y 轴刻度线
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),  # 添加 X 轴刻度线
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
print(tt)
ggsave("../results/TT.png", plot = tt, 
       width = 6, height = 4, dpi = 300)



#E2
results_e2 <- data1 %>%
  group_by(time) %>%
  rstatix::wilcox_test(E ~ group) %>%  
  left_join(
    data1 %>% 
      group_by(time) %>% 
      rstatix::wilcox_effsize(E ~ group),
    by = c("time", "group1", "group2")
  )
print(results_e2)


data1 <- data1 %>% 
  mutate(E = ifelse(E > 4000, NA, E))

summary_data1 <- data1 %>%
  group_by(group, time) %>%
  summarise(
    Mean = mean(E,na.rm = TRUE),
    SD = sd(E,na.rm = TRUE),
    .groups = "drop"
  )



e2<-ggplot() +
  geom_line(
    data = summary_data1,
    aes(x = time, y = Mean, group = group, color = group),
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    data = summary_data1,
    aes(x = time, y = Mean, color = group),
    size = 3,
    position = position_dodge(width = 0.3)
  ) +
  geom_errorbar(
    data = summary_data1,
    aes(x = time, ymin = Mean - SD, ymax = Mean + SD, color = group),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  geom_jitter(
    data = data1,
    aes(x = time, y = E, color = group),
    alpha = 0.4,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.3
    )
  ) +
  labs(
    x = "The diameter of follicles (mm)",
    y = "Serum E2 levels (pg/mL)"
  ) +
  scale_color_manual(
    values = c("expri" = "#E41A1C", "control" = "#377EB8"),
    labels = c("expri" = "mLP", "control" = "GnRH-ant")
  ) +
  #scale_y_continuous(
  # limits = c(0, 6000),
  # breaks = seq(0, 6000, by = 1000),
  # expand = expansion(mult = c(0, 100))
  #) +
  # scale_x_discrete(
  #  labels = c("Initial", "<10", "11-12", "13-15", "16-17", "≥18")
  #) +
  scale_x_continuous(limits = c(0.5,6.5),
                     breaks = c(1,2,3,4,5,6),
                     labels = c("Initial", "<10", "11-12", "13-15", "16-17", "≥18")
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),  # 加粗坐标轴线
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),  # 添加 Y 轴刻度线
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),  # 添加 X 轴刻度线
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
print(e2)
ggsave("../results/E2.png", plot = e2, 
       width = 6, height = 4, dpi = 300)


#LH
results_lh <- data1 %>%
  group_by(time) %>%
  rstatix::wilcox_test(LH ~ group) %>%  # 明确指定包
  left_join(
    data1 %>% 
      group_by(time) %>% 
      rstatix::wilcox_effsize(LH ~ group),
    by = c("time", "group1", "group2")
  )

print(results_lh)

data1 <- data1 %>% 
  mutate(LH = ifelse(LH > 40, NA, LH))


summary_data2 <- data1 %>%
  group_by(group, time) %>%
  summarise(
    Mean = mean(LH,na.rm = TRUE),
    SD = sd(LH,na.rm = TRUE),
    .groups = "drop"
  )


lh<-ggplot() +
  geom_line(
    data = summary_data2,
    aes(x = time, y = Mean, group = group, color = group),
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    data = summary_data2,
    aes(x = time, y = Mean, color = group),
    size = 3,
    position = position_dodge(width = 0.3)
  ) +
  geom_errorbar(
    data = summary_data2,
    aes(x = time, ymin = Mean - SD, ymax = Mean + SD, color = group),
    width = 0.2,
    position = position_dodge(width = 0.3)
  ) +
  geom_jitter(
    data = data1,
    aes(x = time, y = LH, color = group),
    alpha = 0.4,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.3
    )
  ) +
  labs(
    x = "The diameter of follicles (mm)",
    y = "Serum LH levels (IU/L)"
  ) +
  scale_color_manual(
    values = c("expri" = "#E41A1C", "control" = "#377EB8"),
    labels = c("expri" = "mLP", "control" = "GnRH-ant")
  ) +
  scale_y_continuous(
    limits = c(0, 30),
    breaks = seq(0, 30, by = 10),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_discrete(
    labels = c("Initial", "<10", "11-12", "13-15", "16-17", "≥18")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),  # 加粗坐标轴线
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),  # 添加 Y 轴刻度线
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),  # 添加 X 轴刻度线
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
print(lh)
ggsave("../results/LH.png", plot = lh, 
       width = 6, height = 4, dpi = 300)