# Differential Expression Analysis of Dietary Restriction Effects on Mouse Hippocampal Transcriptome

**Course:** BINF6110  
**Author:** Rameen Jafri  
**Dataset:** GSE111778 — Wahl et al. (2018), *Cell Reports*

---

## Background

Caloric restriction (CR) and low-protein, high-carbohydrate (LPHC) diets are among the most well-studied nutritional interventions for extending lifespan and improving healthspan in model organisms. While the systemic metabolic effects of these diets are well characterized, their impact on brain aging at the transcriptomic level remains less understood. This analysis re-examines hippocampal RNA-seq data from Wahl et al. (2018), who compared the effects of caloric restriction and varying levels of dietary protein restriction on gene expression, nutrient-sensing pathways, and cognitive function in 15-month-old mice.

The dataset (GEO accession: GSE111778) contains bulk RNA-seq data from the hippocampus of male and female C57BL/6 mice fed one of five diets from weaning. This analysis includes all five dietary groups: a standard chow control diet (19% protein), a caloric restriction group (ChowCR), and three low-protein, high-carbohydrate diets (15%, 10%, and 5% protein), representing a dose-response gradient of protein restriction. Each group contains 6 biological replicates (3 male, 3 female), for 30 samples total.

---

## Proposed Methods

### 1. Data Acquisition and Quality Control

Raw sequencing reads for all 30 samples (6 per group) were downloaded from the NCBI Sequence Read Archive (SRA accession: SRP135283) using the SRA Toolkit (v3.3.0). Samples correspond to the following dietary groups:

| Group | Description | SRA Accessions |
|-------|-------------|----------------|
| Chow | 19% protein diet | SRR6829578–SRR6829583 |
| ChowCR | Chow + caloric restriction | SRR6829584–SRR6829589 |
| 5% Protein | Low protein diet (5%) | SRR6829590–SRR6829595 |
| 10% Protein | Moderate protein diet (10%) | SRR6829596–SRR6829601 |
| 15% Protein | Moderate protein diet (15%) | SRR6829602–SRR6829607 |

---

#### Total Samples
- 30 RNA-seq samples  
- 5 experimental groups  
- 6 biological replicates per group

Per the recommendations of Williams et al. (2016) and Grigoriev (2020), adapter trimming is not performed prior to pseudoalignment-based quantification, as it has been shown to have minimal impact on downstream gene expression estimates and may reduce mapping rates.

### 2. Quantification

Transcript-level quantification was performed using Salmon (v1.10.3) in selective alignment mode. A transcriptome index was built from the *Mus musculus* GRCm38 (mm10) cDNA reference (Ensembl release 102) using a k-mer length of 31. Each sample was quantified using the following parameters:

```
salmon quant \
    -i reference/salmon_index \
    -l A \
    -r sample.fastq \
    -p 8 \
    --validateMappings \
    -o salmon_output/SAMPLE_ID
```

Transcript-level counts were imported into R using the `tximport` package (v1.36.1), with transcript-to-gene mapping derived from the Ensembl annotation.

Note: while Wahl et al. (2018) used a Tophat/featureCounts pipeline in the original paper, I used Salmon for quantification as it represents current best practice, offering improved speed and accuracy for transcript-level quantification.

### 3. Differential Expression Analysis

Differential expression analysis will be performed using DESeq2 (v1.48.2) in R (v2025.09.2+418). Raw count matrices will be constructed from Salmon output using tximport. Sex will be included as a covariate in the design formula to account for known sex-specific effects on nutrient-sensing pathways reported by Wahl et al. (2018):

```r
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ sex + diet)
```

Each dietary group will be compared against the Chow control:
- Chow vs. ChowCR
- Chow vs. 15% protein
- Chow vs. 10% protein
- Chow vs. 5% protein

Genes with an adjusted p-value (Benjamini-Hochberg correction) < 0.05 and absolute log2 fold change > 1 will be considered differentially expressed.

### 4. Visualization of Data Structure

Overall data structure was assessed using a sample-to-sample distance heatmap computed from variance-stabilizing transformed (VST) counts using the DESeq2 vst() function. Euclidean distances between all 30 samples were calculated and visualized using the pheatmap package, with samples annotated by diet group and sex.

To visualize differential expression patterns, a heatmap of the top 50 most differentially expressed genes from the CR vs Chow comparison (ranked by adjusted p-value) was generated from VST-normalized counts, with rows scaled by z-score to enable cross-gene comparison.

Differential expression results were visualized using volcano plots constructed with ggplot2, displaying log2 fold change against -log10 p-value for each dietary comparison. Genes with adjusted p-value < 0.05 were highlighted, with upregulated genes shown in red and downregulated genes in blue. The top 15 most significant genes were labelled using the ggrepel package.

### 6. Functional Enrichment Analysis – Gene Set Enrichment Analysis (GSEA)

To interpret transcriptional changes at the pathway level, functional annotation and enrichment analysis were performed using Gene Set Enrichment Analysis (GSEA). Ensembl gene identifiers were mapped to Entrez IDs and gene symbols using the org.Mm.eg.db database to enable biological interpretation and compatibility with downstream enrichment tools.

For each dietary comparison (ChowCR, 5% protein, 10% protein, and 15% protein vs Chow), genes were ranked using the Wald test statistic derived from unshrunken DESeq2 results. Genes with missing statistics were excluded, version numbers were removed from Ensembl IDs for consistent mapping, and duplicate Entrez IDs were filtered. The resulting ranked gene list for each contrast was used as input for enrichment analysis.

GSEA was conducted using clusterProfiler::gseGO() with Gene Ontology (GO) Biological Process terms. Gene sets between 15 and 500 genes were considered, and a nominal p-value cutoff of 0.05 was applied. Unlike over-representation analysis, which depends on predefined significance thresholds, GSEA evaluates the entire ranked gene list, allowing detection of coordinated shifts across biologically related pathways.

The following comparisons were analyzed:
- **CR vs Chow**
- **5% Protein vs Chow**
- **10% Protein vs Chow**
- **15% Protein vs Chow**

Significantly enriched pathways were visualized using dotplots, with positive and negative enrichment indicating coordinated upregulation and downregulation, respectively.

---

## Expected Results

Based on the findings of Wahl et al. (2018) and the broader dietary restriction literature, the following outcomes are anticipated:

**Differential expression:** CR and LPHC diets were previously associated with hundreds of differentially expressed genes in the hippocampus compared to control diet. We expect to recapitulate a dose-dependent pattern of gene expression change across the three protein restriction groups, with the 5% protein group showing the greatest divergence from chow controls.

**Pro- and anti-longevity genes:** Based on Wahl et al. (2018), genes in the GenAge database associated with longevity pathways (SIRT1, mTOR, PGC1α) are expected to show differential expression, particularly in the 5% protein group.

**Neurological pathways:** Dendrite morphogenesis, synapse functioning, and neuronal development pathways were previously identified as significantly affected by dietary protein restriction. We expect ORA/GSEA to similarly implicate these GO terms.

**Inflammatory pathways:** Anti-inflammatory cytokine signaling (particularly IL-10 related pathways) may be enriched in the low-protein groups based on prior findings.

**Sex differences:** The original study identified sex-specific effects on nutrient-sensing proteins. If sex is included as a covariate in the DESeq2 model, we may observe reduced noise and improved statistical power.

---

## Results
FILL IN

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| SRA Toolkit | 3.3.0 | Raw data download |
| Salmon | 1.10.3 | Transcript quantification |
| R | v2025.09.2+418 | Statistical analysis |
| DESeq2 | v1.48.2 | Differential expression |
| tximport | v1.36.1 | Count matrix import |
| pheatmap | v1.0.13 | Heatmap visualization |
| ggplot2 | v4.0.2 | General visualization |

---

## References

- Wahl D, et al. (2018). Comparing the Effects of Low-Protein and High-Carbohydrate Diets and Caloric Restriction on Brain Aging in Mice. *Cell Reports*, 25(8), 2234–2243.
- Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
- Patro R, et al. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14, 417–419.
- Soneson C, Love MI, Robinson MD. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521.
- Yu G, et al. (2012). clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS*, 16(5), 284–287.
