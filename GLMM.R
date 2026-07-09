library(raster)
library(vegan)
library(magrittr)
library(data.table)
library(ggplot2)
#library(rgdal)
library(lme4)
library(dplyr)
setwd("D:\\Work\\GO-validation\\test")
#site
site_1900s <- fread("shp/1900s.site",header=T)
site_2000s <- fread("shp/2000s.site",header=T)
site_comb <- rbind(site_1900s,site_2000s)
site_comb <- site_comb[, site := letters[.GRP], by = .(x, y)]
site_comb$site <- as.factor(site_comb$site)
site_comb$ID <- 1:nrow(site_comb)
site_comb$ID <- as.factor(site_comb$ID)
site_comb$Era <- c(rep("1900s",nrow(site_1900s)),rep("2000s",nrow(site_2000s)))
site_comb$Era <- as.factor(site_comb$Era)

#env
env_1900s <- stack("env/bioclim_1900s_ens.tif")
env_2000s <- stack("env/bioclim_2000s_ens.tif")

#extract env
#1900s
extr_1900s <- terra::extract(env_1900s,site_1900s,method="bilinear",ID=F)
extr_1900s_scale <- scale(extr_1900s, center=TRUE, scale=TRUE)
scale_env_1900s <- attr(extr_1900s_scale, 'scaled:scale')
center_env_1900s <- attr(extr_1900s_scale, 'scaled:center')

extr_2000s <- terra::extract(env_2000s,site_2000s,method="bilinear",ID=F)
extr_2000s_scale <- scale(extr_2000s, center=TRUE, scale=TRUE)
scale_env_2000s <- attr(extr_2000s_scale, 'scaled:scale')
center_env_2000s <- attr(extr_2000s_scale, 'scaled:center')

extr_comb <- rbind(extr_1900s,extr_2000s)
extr_comb_scale <- scale(extr_comb, center=TRUE, scale=TRUE)
#load genotype data
rda.1900s.o <- fread("RDA-GO/1900s.overlap.pRDA.genotype")
rda.2000s.o <- fread("RDA-GO/2000s.overlap.pRDA.genotype")

geno_comb <- rbind(rda.1900s.o,rda.2000s.o)
#GLMM测试两个时期的GEA是否相等
df_all <- cbind(site_comb,extr_comb_scale,geno_comb )
df_all$Era <- relevel(df_all$Era, ref = "1900s")

results_df <- data.frame(
  SNP = character(),
  P_value_Interaction = numeric(),
  Slope_1900s = numeric(),
  Slope_2000s = numeric(),
  inter_factor = numeric(),
  stringsAsFactors = FALSE
)

snp_columns <- paste0("V", 1:151)
for (snp in snp_columns) {
  
  # 构建二项分布响应变量矩阵：
  # SNP 值为 0 -> (0, 2)
  # SNP 值为 1 -> (1, 1)
  # SNP 值为 2 -> (2, 0)
  Count_Alt <- df_all[[snp]]
  Count_Ref <- 2 - df_all[[snp]]
  y_matrix <- cbind(Count_Alt, Count_Ref)
  
  # 构建带有随机效应的二项分布混合模型
  model <- glmer(y_matrix ~ bio19 * Era + x + y + (1 | site), #分别对每个bio因子进行处理
                 data = df_all, 
                 family = binomial)
  
  # 显著性 (使用 Anova 提取交互项的 P 值)
  model_anova <- car::Anova(model, type = "III")
  p_val_interaction <- model_anova["bio19:Era", "Pr(>Chisq)"]
  
  # 1900s 的斜率 
  slope_1900s <- fixef(model)["bio19"]
  # 2000s 的斜率 = 1900s斜率 + 交互项系数
  slope_2000s <- slope_1900s + fixef(model)["bio19:Era2000s"]
  
  # 保存结果
  results_df <- rbind(results_df, data.frame(
    SNP = snp,
    P_value_Interaction = p_val_interaction,
    Slope_1900s = slope_1900s,
    Slope_2000s = slope_2000s,
    inter_factor = fixef(model)["bio19:Era2000s"]
  ))
}
#控制假阳性率 (FDR)
results_df$FDR <- p.adjust(results_df$P_value_Interaction, method = "fdr")

# 查看有多少个 SNP 表现出了显著的“空间-时间斜率不一致”
mismatch_snps <- filter(results_df, FDR < 0.05)
print(paste("打破 SFTS 假设的 SNP 数量:", nrow(mismatch_snps)))
#分别保存
bio05.df <- results_df
bio06.df <- results_df
bio18.df <- results_df
bio19.df <- results_df

