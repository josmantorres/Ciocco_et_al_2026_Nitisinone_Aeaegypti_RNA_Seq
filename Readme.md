# RNA-seq analysis of nitisinone treatment in *Aedes aegypti*

## Overview

This repository contains the R scripts and input/output files used for the analysis of an RNA-seq experiment investigating the transcriptional effects of **nitisinone treatment** in *Aedes aegypti* at two time points following treatment.

Nitisinone was administered together with a blood meal. Control insects received a blood meal containing **PBS instead of nitisinone**. The experiment included two independent treatment groups at each time point:

* **Control:** blood meal containing PBS
* **Treated:** blood meal containing nitisinone

RNA was collected at:

* **24 hours after treatment**
* **48 hours after treatment**

The analysis was performed using **DESeq2** and **HTSFilter** in R.

---

## Experimental design

The experiment follows a **2 × 2 factorial design** with two experimental factors:

| Factor                  | Levels           |
| ----------------------- | ---------------- |
| Treatment (`condition`) | Control, Treated |
| Time (`time`)           | 24 h, 48 h       |

The resulting experimental groups are:

| Time | Condition | Treatment  |
| ---- | --------- | ---------- |
| 24 h | Control   | PBS        |
| 24 h | Treated   | Nitisinone |
| 48 h | Control   | PBS        |
| 48 h | Treated   | Nitisinone |

The experimental groups consist of **independent insects**.

---

## Differential expression analysis workflow

The differential expression analysis consists of the following main steps:

1. Import raw gene-level read counts.
2. Import experimental metadata.
3. Construct a DESeq2 dataset.
4. Filter low-expression/high-variability genes using HTSFilter.
5. Generate rlog-transformed expression values.
6. Fit the full factorial DESeq2 model.
7. Test the treatment × time interaction using a likelihood ratio test (LRT).
8. Perform treatment-specific differential expression analyses at 24 h and 48 h.
9. Apply log2 fold-change shrinkage using `apeglm`.

---

## Gene filtering

Low-expression and highly variable genes were filtered using the **HTSFilter** package.

The filtering procedure was run with:

```r
HTSFilter(
    dds_all,
    s.len = 100,
    plot = TRUE
)
```

The resulting filtered dataset was used for the downstream DESeq2 analyses.

---

## Data transformation

For exploratory analyses and visualization, the filtered count data were transformed using the regularized log transformation (`rlog`) implemented in DESeq2:

```r
rlog(
    dds_filtered,
    blind = TRUE,
    fitType = "parametric"
)
```

The resulting rlog-transformed expression matrix is provided as an output file.

The rlog transformation was used for visualization and exploratory analyses and was **not used as input for the differential expression model**. Differential expression was performed using the original integer count data after filtering.

---

# Differential expression analysis

## Full factorial DESeq2 model

Differential expression was initially analyzed using the following model:

```r
~ condition + time + condition:time
```

This model can be expressed as:

```text
Gene expression ~ Treatment + Time + Treatment × Time
```

where:

* `condition` represents nitisinone treatment versus PBS control.
* `time` represents 24 h versus 48 h.
* `condition: time` represents the interaction between treatment and time.

This model allows the independent effects of treatment and time to be separated from their interaction.

---

## Likelihood Ratio Test (LRT)

The treatment × time interaction was evaluated using a likelihood ratio test.

The full model was:

```r
~ condition + time + condition:time
```

and was compared with the reduced model:

```r
~ condition + time
```

The analysis was performed using:

```r
dds_LRT <- DESeq(
    dds_complete,
    test = "LRT",
    reduced = ~ time + condition
)
```

### Interpretation

The LRT tests whether including the `condition:time` interaction significantly improves the explanation of gene expression compared with a model containing only the main effects of treatment and time.

Therefore, genes identified by the LRT represent genes for which the **effect of nitisinone is dependent on time**.

For example, a significant interaction could arise if:

* a gene responds to nitisinone at 24 h but not at 48 h;
* a gene responds at 48 h but not at 24 h;
* the magnitude of the response changes between time points;
* or the direction of the treatment effect changes between time points.

The LRT therefore addresses a different biological question from the individual 24 h and 48 h treatment comparisons.

---

# Treatment-specific comparisons

In addition to the global factorial analysis, treatment effects were evaluated separately at each time point.

## 24 h: Treated vs Control

Samples collected 24 h after treatment were extracted from the complete dataset and analyzed using a model containing only treatment:

```r
~ condition
```

The specific contrast was:

```r
condition_Treated_vs_Control
```

This analysis addresses:

> **What is the transcriptional effect of nitisinone 24 h after administration?**

Log2 fold changes were shrunk using the `apeglm` method.

The analysis used:

```r
lfcShrink(
    dds_24h,
    coef = "condition_Treated_vs_Control",
    res = res_24h,
    type = "apeglm",
    lfcThreshold = 1
)
```

---

## 48 h: Treated vs Control

The same procedure was applied independently to samples collected 48 h after treatment.

The specific contrast was:

```r
condition_Treated_vs_Control
```

This analysis addresses:

> **What is the transcriptional effect of nitisinone 48 h after administration?**

Log2 fold changes were again shrunk using `apeglm`.

---

## Why perform both the LRT and the time-specific comparisons?

The LRT and the individual treatment comparisons answer complementary questions.

### LRT

The LRT asks:

> **Does the effect of nitisinone depend on time?**

### 24 h comparison

The 24 h contrast asks:

> **Which genes differ between nitisinone-treated and control insects at 24 h?**

### 48 h comparison

The 48 h contrast asks:

> **Which genes differ between nitisinone-treated and control insects at 48 h?**

Consequently, a gene can show a significant treatment effect at one time point without necessarily showing a significant treatment × time interaction. Conversely, a significant interaction indicates that the treatment effect differs between time points.

---

# Repository structure

A recommended repository structure is:

```text
.
├── README.md
│
├── scripts/
│   └──Read_processing_mapping_counting_details.txt
│   └──RNAseq_DESeq2_nitisinone.R
│
├── input/
│   ├── Counts_all.txt
│   └── Experiment_description_all.txt
│
└── results/
    ├── RLOG_transformed_values_all.txt
    ├── DESeq2_LRT.txt
    ├── DESeq2_statistics_24h.txt
    └── DESeq2_statistics_48h.txt
```

---

# Input files

## `Counts_all.txt`

Raw gene-level RNA-seq read counts.

The first column contains gene identifiers, while the remaining columns correspond to individual RNA-seq samples.

The count matrix is used as input for the DESeq2 analysis.

---

## `Experiment_description_all.txt`

Experimental metadata associated with each RNA-seq sample.

The metadata includes at least:

* `sample`: sample identifier
* `condition`: treatment condition (`Control` or `Treated`)
* `time`: sampling time point (`24` or `48`)

The sample identifiers in the metadata must correspond to the column names of the count matrix.

---

# Output files

## `RLOG_transformed_values_all.txt`

Matrix containing rlog-transformed expression values for genes retained after HTSFilter.

These values are intended for exploratory analyses and visualization.

---

## `DESeq2_LRT.txt`

Results of the likelihood ratio test evaluating the `condition × time` interaction.

This file contains the complete set of genes tested by the LRT.

The adjusted p-value (`padj`) can be used to identify genes for which the transcriptional effect of nitisinone significantly differs between 24 h and 48 h.

---

## `DESeq2_statistics_24h.txt`

Differential expression results for:

```text
24 h: Nitisinone-treated vs PBS control
```

The reported log2 fold changes have been shrunk using `apeglm`.

---

## `DESeq2_statistics_48h.txt`

Differential expression results for:

```text
48 h: Nitisinone-treated vs PBS control
```

The reported log2 fold changes have been shrunk using `apeglm`.

---

# Software and packages

The analysis was performed in R using the following packages:

* **DESeq2** — differential expression analysis of RNA-seq count data.
* **HTSFilter** — filtering of low-expression/high-variability genes.
* **apeglm** — adaptive shrinkage of log2 fold changes.

The exact package versions used for the analysis should be reported here to facilitate reproducibility.

For example:

```text
Rstudio version: 2026.07.1
DESeq2 version: 1.42.1
HTSFilter version: 1.42
apeglm version: 1.24
```

---

# Reproducibility

The complete analysis is implemented in:

```text
scripts/RNAseq_DESeq2_nitisinone.R
```

Before running the script, update the `input_dir` and `output_dir` variables to point to the appropriate directories.

The analysis requires the input count matrix and corresponding experimental metadata described above.

---

# Notes

* Raw read counts, rather than normalized expression values, are used as input for DESeq2.
* The rlog-transformed data are generated for exploratory analyses and visualization only.
* The full factorial model is used to evaluate treatment, time, and treatment × time effects.
* The LRT specifically evaluates the treatment × time interaction.
* Separate models are used to estimate the nitisinone treatment effect at 24 h and 48 h.
* Log2 fold changes from the time-specific comparisons are shrunk using `apeglm`.
* The experimental groups consist of independent insects.

---

# Citation

If you use these analyses or code, please cite the associated publication:

> **[Add manuscript citation here]**

Please also cite the software packages used in the analysis, including DESeq2, HTSFilter, and apeglm.

