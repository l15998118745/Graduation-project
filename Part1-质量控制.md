# SLE cfRNA 数据分析工作流程（当前版：TPM 数据）

> 项目：毕业设计（cfRNA / SLE 方向）  
> 数据：3069 个样本 × 19665 个 mRNA 基因的表达矩阵 + 3069 行 metadata  
> 当前依据：`Part1-质量控制与PCA.R`  
> 标记说明：正常文字为当前流程；~~删除线~~为旧流程中保留备查、但现在不再执行的步骤。

---

## 总览

```text
TPM 表达矩阵 + metadata
  │
  ├─ 1. 环境准备与读入数据
  ├─ 2. 样本 ID 和行列名匹配
  ├─ 3. 数值及数据格式诊断
  ├─ 4. 转为“基因 × 样本”矩阵
  ├─ 5. log2(TPM + 1)
  ├─ 6. 清理 metadata、提取批次
  ├─ 7. 初始 PCA
  ├─ 8. adonis2 + MANOVA 定量数据结构
  ├─ 9. 线性模型去除批次效应
  └─ 10. 批次校正后 PCA 与 MANOVA
```

旧版归一化分支已经停用：

~~expected counts → 取整 → DESeq2 对象 → 基因预过滤 → size factor 归一化 → VST → VST-PCA~~

**当前原则**：输入矩阵已经完成 TPM 归一化；不再使用 DESeq2 重新估计 size factor，也不做 VST。

---

## 1. 环境准备与读入数据

```r
library(data.table)
library(ggplot2)
library(matrixStats)
library(irlba)
library(vegan)

dir.create("results", showWarnings = FALSE)

mx   <- fread("data/expression_matrix aligned_mRNA_matrix")
meta <- fread("data/metadata anno_final_aligned")

cat("mx 维度:", paste(dim(mx), collapse = " x "), "\n")
cat("meta 维度:", paste(dim(meta), collapse = " x "), "\n")
```

| 对象 | 预期维度 | 内容 |
| --- | ---: | --- |
| `mx` | 3069 × 19666 | `sample_id` + 19665 个基因 |
| `meta` | 3069 × 27 | 样本注释 |

> ~~旧流程需要加载 `DESeq2` 并构建 DESeq2 对象。~~ 当前 Part1 不再需要。

---

## 2. 样本 ID 与行列名匹配

**目的**：确保表达矩阵与 metadata 一一对应，避免 PCA 分数和样本标签错配。

```r
stopifnot(
  all(mx$sample_id %in% meta$id),
  nrow(meta) == nrow(mx),
  !anyDuplicated(mx$sample_id),
  !anyDuplicated(meta$id)
)
```

当前定案：`mx$sample_id` 对应 `meta$id`；`meta$sample` 不是本流程的文库级匹配键。

---

## 3. 数值与数据格式诊断

```r
x <- as.matrix(mx[, -1])
storage.mode(x) <- "double"

min(x, na.rm = TRUE)
max(x, na.rm = TRUE)
mean(is.na(x))
mean(x == 0, na.rm = TRUE)

rs <- rowSums(x, na.rm = TRUE)
summary(rs)

gm <- colMeans(x, na.rm = TRUE)
gv <- colVars(x, na.rm = TRUE)
ok <- gm > 0 & gv > 0 & is.finite(gm) & is.finite(gv)
fit <- lm(log10(gv[ok]) ~ log10(gm[ok]))
coef(fit)[2]
```

均值—方差图保存为 `results/mean_var.png`。

### 当前数据类型结论

矩阵是 **mRNA 与其他 RNA 种类一起完成 TPM 后，再提取出的 mRNA 子集**。所以：

- 数值按 TPM 表达量处理；
- 每个样本只保留 mRNA，行和不必等于 1,000,000；
- 行和不等于 1,000,000 不能据此把矩阵重新判定为原始 count；
- 当前 PCA 使用 `log2(TPM + 1)`，不再做额外的文库大小归一化。

~~旧笔记把矩阵判定为 expected counts，并认为可以直接进入 DESeq2。~~  
~~当前不再沿用该判断，也不把 TPM 值取整后当作 count。~~

---

## 4. 整理表达矩阵并进行对数转换

```r
expr <- t(as.matrix(mx[, -1]))
colnames(expr) <- mx$sample_id
rownames(expr) <- colnames(mx)[-1]
expr[is.na(expr)] <- 0

logmat <- log2(expr + 1)
rm(x, expr)
```

`logmat` 为 19665 个基因 × 3069 个样本，是初始 PCA、批次校正和校正后 PCA 的工作矩阵。

### 已停用的旧预处理

~~将表达值 `round()` 取整并强制转换为 integer。~~  
~~使用 `DESeqDataSetFromMatrix()` 构建 `design = ~1` 的 DESeq2 对象。~~  
~~保留“至少在 300 个样本中 count ≥ 5”的基因，将 19665 个基因过滤为 4374 个。~~

这些步骤都依赖“输入是 expected counts”的前提，与当前 TPM 流程不一致。

---

## 5. 归一化步骤变更记录

### 当前做法

输入数据已经完成 TPM 归一化。当前只进行：

```r
logmat <- log2(expr + 1)
```

`log2(TPM + 1)` 是对数转换，用于压缩高表达值尺度、改善 PCA 展示；它不是重新归一化。

### 已停用做法

~~`estimateSizeFactors(dds_all, type = "poscounts")`：估计样本 size factor。~~  
~~`varianceStabilizingTransformation(dds_all, blind = TRUE)`：进行 VST。~~  
~~保存 `results/mat_vst.rds`，以 `mat_vst` 作为 PCA 输入。~~  
~~基于 VST 后的行方差选择 top 5000 高变基因。~~

停用原因：这些操作属于旧版 expected-count/DESeq2 流程；当前矩阵已是 TPM 数据。

---

## 6. 清理 metadata 并提取批次

```r
meta[, cell_type_clean := toupper(cell_type)]

get_batch <- function(x) {
  out <- character(length(x))
  has_date <- grepl("^\\d{8}", x)
  out[has_date] <- substr(x[has_date], 1, 8)
  out[!has_date] <- sub("_.*$", "", x[!has_date])
  out
}
```

- `toupper()` 统一 `Neu/NEU`、`Plt/PLT` 等大小写差异。
- ID 以 8 位日期开头时，批次取该日期。
- 没有日期前缀时，批次取第一个下划线前的字段。

---

## 7. 初始 PCA：观察未经批次校正的数据结构

```r
pca <- prcomp_irlba(t(logmat), n = 10, center = TRUE)
var_exp0 <- round(100 * pca$sdev[1:5]^2 / pca$totalvar, 1)

pc_df <- as.data.frame(pca$x[, 1:3])
pc_df$sample <- colnames(logmat)

meta_key <- meta[, .(id, cell_type_clean, BinaryGroup, patient)]
setnames(meta_key, "id", "sample")
pc_df <- merge(pc_df, meta_key, by = "sample")
saveRDS(pc_df, "results/pc_df.rds")
```

合并后应仍为 3069 行。当前脚本输出：

- `results/pca_by_celltype.pdf`
- `results/pca_by_group.pdf`

脚本随后重新提取前 10 个 PC，用于结构检验，并输出：

- `results/pc_df2.rds`
- `results/pca_log_celltype.pdf`
- `results/pca_log_group.pdf`
- `results/pca_log_batch.pdf`

---

## 8. 定量判断数据结构的主导因素

### adonis2 置换检验

```r
d <- dist(pc_df[, paste0("PC", 1:10)])
adonis2(d ~ cell_type_clean, data = pc_df, permutations = 199)
adonis2(d ~ batch,            data = pc_df, permutations = 199)
adonis2(d ~ BinaryGroup,      data = pc_df, permutations = 199)
```

### MANOVA

```r
Y <- as.matrix(pc_df[, paste0("PC", 1:10)])
summary(manova(Y ~ cell_type_clean, data = pc_df), test = "Wilks")
summary(manova(Y ~ batch,            data = pc_df), test = "Wilks")
summary(manova(Y ~ BinaryGroup,      data = pc_df), test = "Wilks")
```

当前脚本同时保留置换检验和 MANOVA，用来比较细胞类型、批次和疾病分组与前 10 个 PC 整体结构的关系。

---

## 9. 线性模型批次校正

**目的**：去除批次对应的线性效应，同时保护细胞类型与疾病分组信号。

```r
m_mat <- match(colnames(logmat), meta$id)
meta_mat <- meta[m_mat, ]
meta_mat$batch <- get_batch(meta_mat$id)

X <- model.matrix(~ batch + cell_type_clean + BinaryGroup, data = meta_mat)
colnames(X) <- make.names(colnames(X))

beta <- solve(crossprod(X), crossprod(X, t(logmat)))
b_coef <- beta[grep("^batch", rownames(beta)), , drop = FALSE]

mat_cor <- logmat - t(
  X[, grep("^batch", colnames(X)), drop = FALSE] %*% b_coef
)
```

模型同时纳入 batch、细胞类型和疾病分组，但只减去 batch 的拟合部分。这是**批次效应校正**，不是 TPM 或 size-factor 归一化。

---

## 10. 批次校正后 PCA 与复检

```r
pca3 <- prcomp_irlba(t(mat_cor), n = 10, center = TRUE)
pc_df3 <- as.data.frame(pca3$x[, 1:10])
pc_df3$sample <- colnames(mat_cor)
pc_df3 <- merge(pc_df3, meta_key, by = "sample")

var_exp3 <- 100 * pca3$sdev[1:10]^2 / pca3$totalvar
```

输出：

- `results/pca_log_celltype_corrected.pdf`
- `results/pca_log_group_corrected.pdf`

校正后复检：

```r
Y3 <- as.matrix(pc_df3[, paste0("PC", 1:10)])
summary(manova(Y3 ~ cell_type_clean, data = pc_df3), test = "Wilks")
summary(manova(Y3 ~ BinaryGroup,      data = pc_df3), test = "Wilks")
```

---

## 关键文件清单

| 文件 | 内容 | 状态 |
| --- | --- | --- |
| `data/expression_matrix aligned_mRNA_matrix` | TPM 表达矩阵 | 当前输入 |
| `data/metadata anno_final_aligned` | 样本注释 | 当前输入 |
| `results/mean_var.png` | 均值—方差诊断图 | 当前输出 |
| `results/pc_df.rds` | 初始 PCA 前 3 个 PC + metadata | 当前输出 |
| `results/pc_df2.rds` | 初始 PCA 前 10 个 PC + metadata + batch | 当前输出 |
| `results/pca_log_*.pdf` | 初始 log2-TPM PCA 图 | 当前输出 |
| `results/pca_log_*_corrected.pdf` | 批次校正后 PCA 图 | 当前输出 |
| ~~`results/mat_vst.rds`~~ | ~~VST 表达矩阵~~ | ~~旧流程，已停用~~ |
| ~~`results/pca2_vst.rds`~~ | ~~VST-PCA 对象~~ | ~~旧流程，已停用~~ |
| ~~`results/pca_vst_*.pdf`~~ | ~~VST-PCA 图~~ | ~~旧流程，已停用~~ |

---

## 当前已定案的技术结论

1. 数据按 **TPM** 处理；它是所有 RNA 一起归一化后提取的 mRNA 子集，因此 mRNA 行和不要求等于 1,000,000。
2. 合并键统一使用 `mx$sample_id` 对应 `meta$id`。
3. 当前工作矩阵是 `log2(TPM + 1)`，不再做 DESeq2 size-factor 归一化或 VST。
4. 初始 PCA 用于观察未经校正的整体结构。
5. `adonis2` 和 MANOVA 用于比较细胞类型、批次和疾病分组的影响。
6. 批次校正采用线性模型，并保护 `cell_type_clean` 和 `BinaryGroup`。
7. 校正后重新做 PCA 和 MANOVA，检查生物学结构是否保留。

## 后续分析提醒

- 当前 Part1 解决质量检查、结构探索和 PCA 层面的批次校正。
- 是否可以直接用 TPM 做差异表达取决于后续统计方法；不要把 TPM 取整后直接当原始 count 输入 DESeq2。
- 正式差异分析仍需考虑同一受试者的重复测量、细胞类型分层，以及 batch 与生物学分组是否混杂。
