## =====================================================================
## Part1  质量控制与 PCA(TPM 数据版)
## 背景: 该矩阵是 mRNA 与其他 RNA 种类一起做完 TPM 后抽出的 mRNA 子集,
##       所以每行(样本)总和不是 1e6, 而是该样本 mRNA 占总 RNA 的份额。
## 思路: 读入 → 格式/数值诊断 → log2(TPM+1) PCA 看整体结构
##       → 置换检验找主导因素 → 线性批次校正后再看 PCA
## 运行: 工作目录设为项目根目录, 然后
##       Rscript scripts/Part1-质量控制与PCA.R
## 依赖: data.table, ggplot2, matrixStats, irlba, vegan
## =====================================================================

library(data.table)
library(ggplot2)
library(matrixStats)
library(irlba)

dir.create("results", showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1. 读入数据
## ---------------------------------------------------------------------
mx   <- fread("data/expression_matrix aligned_mRNA_matrix")  # 3069 样本 x 19666 列(sample_id + 19665 基因)
meta <- fread("data/metadata anno_final_aligned")            # 3069 x 27

cat("mx 维度:", paste(dim(mx), collapse = " x "), "\n")
cat("meta 维度:", paste(dim(meta), collapse = " x "), "\n")

## 矩阵 sample_id 与 metadata 的 id 列一一对应(注意不是 sample 列)
stopifnot(
  all(mx$sample_id %in% meta$id),
  nrow(meta) == nrow(mx),
  !anyDuplicated(mx$sample_id), !anyDuplicated(meta$id)
)
cat("样本ID全部匹配 metadata$id:", nrow(mx), "/", nrow(mx), "\n")

## ---------------------------------------------------------------------
## 2. 数值诊断
## ---------------------------------------------------------------------
x <- as.matrix(mx[, -1])          # 3069 x 19665, 约 480MB
storage.mode(x) <- "double"

cat("\n--- 数值特征 ---\n")
cat("最小值:", min(x, na.rm = TRUE), " 最大值:", max(x, na.rm = TRUE), "\n")
cat("缺失占比:", mean(is.na(x)), "\n")
cat("零值占比(该基因在该样本无读段):", mean(x == 0, na.rm = TRUE), "\n")

rs <- rowSums(x, na.rm = TRUE)    # 每样本 mRNA TPM 总量(占总 RNA 的份额)
cat("\n--- 每样本 mRNA TPM 总量(rowSums) ---\n"); print(summary(rs))
cat("(若为单独归一化的 TPM, 这里应恒为 1e6;\n")
cat(" 现在各行不等, 正是'与所有 RNA 一起归一化'的特征)\n")

## 均值-方差关系(对数刻度看趋势)
gm <- colMeans(x, na.rm = TRUE)
gv <- colVars(x, na.rm = TRUE)
ok <- gm > 0 & gv > 0 & is.finite(gm) & is.finite(gv)
fit <- lm(log10(gv[ok]) ~ log10(gm[ok]))
cat("\n均值-方差 log-log 拟合斜率:", round(coef(fit)[2], 2), "\n")

png("results/mean_var.png", width = 800, height = 600)
plot(log10(gm[ok]), log10(gv[ok]), pch = 16, cex = 0.3,
     xlab = "log10 mean", ylab = "log10 var",
     main = paste0("slope = ", round(coef(fit)[2], 2)))
abline(fit, col = "red", lwd = 2)
dev.off()

## ---------------------------------------------------------------------
## 3. 整理表达矩阵(基因 x 样本)并取 log2(TPM+1)
## ---------------------------------------------------------------------
expr <- t(as.matrix(mx[, -1]))
colnames(expr) <- mx$sample_id
rownames(expr) <- colnames(mx)[-1]
expr[is.na(expr)] <- 0
cat("\nexpr 维度:", paste(dim(expr), collapse = " x "), "\n")

logmat <- log2(expr + 1)          # TPM 的工作刻度, 之后的 PCA/limma 都用它
rm(x, expr)

## ---------------------------------------------------------------------
## 4. 清理 metadata 标签
## ---------------------------------------------------------------------
meta[, cell_type_clean := toupper(cell_type)]  # Neu->NEU, Plt->PLT
print(sort(table(meta$cell_type_clean)))
print(table(meta$BinaryGroup))
cat("受试者总数:", length(unique(meta$patient)), "\n")

## 批次 = 样本ID的 8 位测序日期前缀(如 20230606); 没有的取第一个下划线前段
get_batch <- function(x) {
  out <- character(length(x))
  has_date <- grepl("^\\d{8}", x)
  out[has_date] <- substr(x[has_date], 1, 8)
  out[!has_date] <- sub("_.*$", "", x[!has_date])
  out
}

## ---------------------------------------------------------------------
## 5. log2(TPM+1) PCA(不做其他处理的初版, 看初始结构)
## ---------------------------------------------------------------------
pca <- prcomp_irlba(t(logmat), n = 10, center = TRUE)   # irlba 只要前10个PC
var_exp0 <- round(100 * pca$sdev[1:5]^2 / pca$totalvar, 1)
cat("\nlog2-PCA 前5个主成分解释方差(%):", var_exp0, "\n")

pc_df <- as.data.frame(pca$x[, 1:3])
pc_df$sample <- colnames(logmat)   # irlba 不打乱行序, 直接对应
meta_key <- meta[, .(id, cell_type_clean, BinaryGroup, patient)]
setnames(meta_key, "id", "sample")
pc_df <- merge(pc_df, meta_key, by = "sample")
cat("合并后 pc_df 行数:", nrow(pc_df), "(应为 3069)\n")
saveRDS(pc_df, "results/pc_df.rds")

p1 <- ggplot(pc_df, aes(PC1, PC2, color = cell_type_clean)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "PCA colored by cell type (log2 TPM)") + theme_bw()
p2 <- ggplot(pc_df, aes(PC1, PC2, color = BinaryGroup)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = "PCA colored by group (log2 TPM)") + theme_bw()
ggsave("results/pca_by_celltype.pdf", p1, width = 9, height = 6)
ggsave("results/pca_by_group.pdf",    p2, width = 8, height = 6)
rm(pca)

## ---------------------------------------------------------------------
## 6. 哪一因素主导数据结构: adonis2 置换检验 + 不置换的 MANOVA
## ---------------------------------------------------------------------
library(vegan)

pc_df$batch <- get_batch(pc_df$sample)
npc <- 10
pc_full <- as.data.frame(pca_full <- NULL)   # 占位防误读
## 需要前10个PC: 重新取(上面只留了3列)
pca10 <- prcomp_irlba(t(logmat), n = 10, center = TRUE)
pc_df <- as.data.frame(pca10$x[, 1:npc])
pc_df$sample <- colnames(logmat)
pc_df <- merge(pc_df, meta_key, by = "sample")
pc_df$batch <- get_batch(pc_df$sample)
saveRDS(pc_df, "results/pc_df2.rds")

var_exp <- 100 * pca10$sdev[1:npc]^2 / pca10$totalvar
pc_lab  <- paste0("PC", 1:npc, " (", round(var_exp, 1), "%)")
cat("前5个PC解释方差(%):", round(var_exp[1:5], 1), "\n")

p1 <- ggplot(pc_df, aes(PC1, PC2, color = cell_type_clean)) +
  geom_point(size = 1, alpha = 0.6) + theme_bw() +
  labs(x = pc_lab[1], y = pc_lab[2], title = "log2-PCA by cell type")
p2 <- ggplot(pc_df, aes(PC1, PC2, color = BinaryGroup)) +
  geom_point(size = 1, alpha = 0.6) + theme_bw() +
  labs(x = pc_lab[1], y = pc_lab[2], title = "log2-PCA by group")
p3 <- ggplot(pc_df, aes(PC1, PC2, color = batch)) +
  geom_point(size = 1, alpha = 0.6) + theme_bw() +
  labs(x = pc_lab[1], y = pc_lab[2], title = "log2-PCA by batch")
ggsave("results/pca_log_celltype.pdf", p1, width = 9, height = 6)
ggsave("results/pca_log_group.pdf",    p2, width = 8, height = 6)
ggsave("results/pca_log_batch.pdf",    p3, width = 9, height = 6)

d <- dist(pc_df[, paste0("PC", 1:npc)])
print(adonis2(d ~ cell_type_clean, data = pc_df, permutations = 199))
print(adonis2(d ~ batch,            data = pc_df, permutations = 199))
print(adonis2(d ~ BinaryGroup,      data = pc_df, permutations = 199))

Y <- as.matrix(pc_df[, paste0("PC", 1:npc)])
print(summary(manova(Y ~ cell_type_clean, data = pc_df), test = "Wilks"))
print(summary(manova(Y ~ batch,            data = pc_df), test = "Wilks"))
print(summary(manova(Y ~ BinaryGroup,      data = pc_df), test = "Wilks"))

## ---------------------------------------------------------------------
## 7. 批次校正后 PCA(线性模型校正, 保护细胞类型与分组变量)
## ---------------------------------------------------------------------
m_mat <- match(colnames(logmat), meta$id)
meta_mat <- meta[m_mat, ]
meta_mat$batch <- get_batch(meta_mat$id)
X <- model.matrix(~ batch + cell_type_clean + BinaryGroup, data = meta_mat)
colnames(X) <- make.names(colnames(X))
beta <- solve(crossprod(X), crossprod(X, t(logmat)))   # 每个基因回归
b_coef <- beta[grep("^batch", rownames(beta)), , drop = FALSE]
mat_cor <- logmat - t(X[, grep("^batch", colnames(X)), drop = FALSE] %*% b_coef)

pca3 <- prcomp_irlba(t(mat_cor), n = 10, center = TRUE)
pc_df3 <- as.data.frame(pca3$x[, 1:10])
pc_df3$sample <- colnames(mat_cor)
pc_df3 <- merge(pc_df3, meta_key, by = "sample")

var_exp3 <- 100 * pca3$sdev[1:10]^2 / pca3$totalvar
pc_lab3  <- paste0("PC", 1:10, " (", round(var_exp3, 1), "%)")

p1c <- ggplot(pc_df3, aes(PC1, PC2, color = cell_type_clean)) +
  geom_point(size = 1, alpha = 0.6) + theme_bw() +
  labs(x = pc_lab3[1], y = pc_lab3[2],
       title = "Batch-corrected PCA by cell type")
p2c <- ggplot(pc_df3, aes(PC1, PC2, color = BinaryGroup)) +
  geom_point(size = 1, alpha = 0.6) + theme_bw() +
  labs(x = pc_lab3[1], y = pc_lab3[2],
       title = "Batch-corrected PCA by group")
ggsave("results/pca_log_celltype_corrected.pdf", p1c, width = 9, height = 6)
ggsave("results/pca_log_group_corrected.pdf",    p2c, width = 8, height = 6)

## 校正后各因素还显不显著(快速 MANOVA)
Y3 <- as.matrix(pc_df3[, paste0("PC", 1:10)])
print(summary(manova(Y3 ~ cell_type_clean, data = pc_df3), test = "Wilks"))
print(summary(manova(Y3 ~ BinaryGroup,      data = pc_df3), test = "Wilks"))

cat("\nPart1 完成, 结果在 results/ 下\n")
