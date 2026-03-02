# Differential Expression Analysis of Dietary Restriction Effects on Mouse Hippocampal Transcriptome

**Course:** BINF6110  
**Author:** Rameen Jafri  
**Dataset:** GSE111778 — Wahl et al. (2018), *Cell Reports*

---

## Background

Aging is the primary risk factor for a broad range of neurodegenerative and cognitive disorders, and identifying interventions that slow brain aging at the molecular level is a major goal of biomedical research. Dietary restriction is one of the most robust and evolutionarily conserved interventions known to extend lifespan and improve healthspan across model organisms ranging from yeast to primates (Fontana & Partridge, 2015). Two of the most well-studied dietary interventions are caloric restriction (CR), which reduces total energy intake without malnutrition, and low-protein, high-carbohydrate (LPHC) diets, which manipulate macronutrient composition rather than total caloric intake. While both interventions have been shown to extend lifespan and modulate nutrient-sensing pathways including mTOR, SIRT1, and IGF-1, their effects on brain aging at the transcriptomic level remain poorly characterised.

The hippocampus is a particularly important brain region in the context of aging, as it is critical for memory consolidation, spatial navigation, and cognitive function, and shows early and pronounced vulnerability to age-related decline. Understanding how dietary interventions modulate hippocampal gene expression may therefore provide insight into the molecular mechanisms underlying diet-induced improvements in cognitive aging. Wahl et al. (2018) addressed this question by profiling hippocampal gene expression in 15-month-old mice fed one of five diets from weaning, finding that both CR and LPHC diets altered expression of genes involved in nutrient-sensing, circadian rhythm, and neuronal function. However, their analysis used a TopHat2/featureCounts alignment pipeline with FPKM normalisation, which has since been superseded by more accurate pseudoalignment-based approaches.

This analysis re-examines the same dataset (GEO accession: GSE111778) using a modern RNA-seq pipeline (Salmon pseudoalignment with DESeq2 differential expression analysis) to assess whether the original findings are reproducible with current best-practice methods, and to identify additional pathway-level insights through Gene Set Enrichment Analysis (GSEA). The dataset contains bulk RNA-seq data from the hippocampus of male and female C57BL/6 mice across five dietary groups: a standard chow control diet (19% protein), a caloric restriction group (ChowCR), and three LPHC diets (15%, 10%, and 5% protein) representing a dose-response gradient of protein restriction. Each group contains 6 biological replicates (3 male, 3 female), for 30 samples total. Sex was included as a covariate in all analyses to account for known sex-specific effects on nutrient-sensing pathways reported by Wahl et al. (2018).

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

Note: Salmon was selected over alignment-based approaches such as TopHat2 (used by Wahl et al.) because pseudoalignment-based quantification has been shown to produce comparable or superior accuracy at a fraction of the computational cost (Patro et al., 2017). Unlike traditional aligners such as TopHat2 or STAR, which map reads to the entire genome including introns and intergenic regions, Salmon maps directly to the known transcriptome, making it faster and less memory-intensive while avoiding spurious alignments to non-transcribed regions. Transcript-level estimates were summarised to gene level using tximport, which has been shown to improve gene-level inference compared to methods that quantify at the gene level directly (Soneson et al., 2015). DESeq2 was preferred over FPKM-based normalisation used in the original study, as FPKM is not suitable for between-sample comparisons and does not model count-based overdispersion, whereas DESeq2's negative binomial model provides better statistical control and more accurate fold change estimation (Love et al., 2014). GSEA was chosen over ORA because it evaluates the entire ranked gene list rather than relying on an arbitrary significance threshold, improving sensitivity for detecting coordinated pathway-level changes (Yu et al., 2012).

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

### Sample Quality and Overall Transcriptional Structure
Sample-to-sample distance heatmap computed from variance-stabilizing transformed (VST) counts across all 30 hippocampal samples (Figure 1) was used to asses sample quality and transcriptional structure. One sample (SRR6829578, male Chow) was immediately identifiable as a clear outlier, exhibiting markedly greater transcriptional distance from all other samples, visible as a distinct bright row and column in the heatmap. Despite this outlier, the sample was retained in the analysis as it passed all technical quality metrics including mapping rate and library size, and its removal was not required by the original study's analytical framework. Beyond the outlier, hierarchical clustering revealed that samples did not cluster primarily by diet group, but instead showed substantial inter-individual variability consistent with the high biological heterogeneity expected in hippocampal tissue. Sex was observed as a notable source of transcriptional variation, with male and female samples showing partial separation in the clustering pattern. Interestingly, Wahl et al. (2018) reported no sex-based separation in their PCA of gene expression data and therefore combined sexes for their transcriptional analysis. However, as the authors did identify significant sex-specific differences in hippocampal protein expression — particularly for SIRT1, mTOR, and PGC1α — and given the clear sex-associated variation visible in our sample distance heatmap, sex was retained as a covariate in the DESeq2 design formula to account for this source of variance.

<img width="2400" height="1687" alt="image" src="https://github.com/user-attachments/assets/7748458c-cdcb-42bc-b363-ca1a45516de5" />

&nbsp;
**Figure 1. Sample-to-sample distance heatmap of hippocampal gene expression across all 30 samples.** 
Euclidean distances were calculated from variance-stabilizing transformed (VST) counts and visualised using hierarchical clustering. Darker blue indicates greater transcriptional similarity between samples; lighter blue indicates greater dissimilarity. Samples are annotated by diet group (Chow, ChowCR, P5, P10, P15) and sex (male, female). 

### Differential Gene Expression 
Differential expression analysis was performed using DESeq2 with sex included as a covariate in the design formula (~ sex + diet), with each dietary group compared against the Chow control diet (19% protein). After low-count filtering — retaining only genes with a minimum of 10 counts in at least 6 samples — 15,133 genes were retained for statistical testing across all comparisons. The number of differentially expressed genes (DEGs, adjusted p-value < 0.05, Benjamini-Hochberg correction) varied substantially across dietary comparisons (Table 1). Caloric restriction produced the strongest transcriptional response with 366 DEGs, followed by 10% protein (116 DEGs), 15% protein (52 DEGs), and 5% protein (4 DEGs).

**Table 1. Summary of differentially expressed genes for each dietary comparison against Chow control.** 

DEGs were identified using DESeq2 with sex as a covariate (design: ~ sex + diet). Significance threshold: adjusted p-value < 0.05 (Benjamini-Hochberg correction). Low-count filtering retained 15,133 genes for testing. CR = caloric restriction; P5, P10, P15 = 5%, 10%, and 15% protein diets respectively.
| Comparison | Total DEGs | Upregulated | Downregulated |
|------------|-----------|-------------|---------------|
| CR vs Chow | 366 | 194 | 172 |
| 5% Protein vs Chow | 4 | 2 | 2 |
| 10% Protein vs Chow | 116 | 70 | 46 |
| 15% Protein vs Chow | 52 | 38 | 14 |

### Top Differentially Expressed Genes
The top 10 most significantly differentially expressed genes for each dietary comparison are summarised in Table 2. Across all three protein restriction groups (5%, 10%, and 15% protein), Gpr17 was consistently the most significantly upregulated gene. In the CR comparison, the top hits included Bmal1 (log2FC = -0.40, padj = 5.84×10⁻⁶), Clock (log2FC = -0.25, padj = 6.93×10⁻⁶), Cited2 (log2FC = 0.48, padj = 6.93×10⁻⁶), and Dbp (log2FC = 0.58, padj = 1.27×10⁻⁵). In the 10% protein comparison, top hits included Gpr17 (log2FC = 0.51, padj = 3.17×10⁻⁶), Banp (log2FC = 0.82, padj = 1.14×10⁻⁴), and Sema4b (log2FC = 0.38, padj = 6.85×10⁻⁴). In the 15% protein comparison, Gpr17 (log2FC = 0.54, padj = 3.05×10⁻⁷) was again the top hit, followed by Plin4 (log2FC = 2.00, padj = 1.75×10⁻⁵) and Hspa5 (log2FC = 0.36, padj = 2.58×10⁻³). These findings closely replicate the top gene results reported by Wahl et al. (2018), with Gpr17, Dbp, Cited2, Sema4b, Hspa5, and Plin4 all appearing in both analyses.

**Table 2. Top 10 differentially expressed genes for each dietary comparison against Chow control.** 

Genes are ranked by adjusted p-value. Log2 fold change values are apeglm-shrunken estimates. padj = Benjamini-Hochberg adjusted p-value.
| Rank | CR Gene | CR log2FC | CR padj | P5 Gene | P5 log2FC | P5 padj | P10 Gene | P10 log2FC | P10 padj | P15 Gene | P15 log2FC | P15 padj |
|------|---------|-----------|---------|---------|-----------|---------|----------|------------|----------|----------|------------|----------|
| 1 | Bmal1 | -0.398 | 5.84e-06 | Gpr17 | 0.396 | 2.74e-03 | Gpr17 | 0.509 | 3.17e-06 | Gpr17 | 0.54 | 3.05e-07 |
| 2 | Clock | -0.245 | 6.93e-06 | Mrpl23-ps1 | 0 | 2.74e-03 | Banp | 0.82 | 1.14e-04 | Plin4 | 1.998 | 1.75e-05 |
| 3 | Cited2 | 0.478 | 6.93e-06 | B3galnt2 | -0.392 | 1.54e-02 | Celf6 | 0.485 | 2.05e-04 | Sdf2l1 | 0.649 | 1.19e-04 |
| 4 | Dbp | 0.582 | 1.27e-05 | Gpr21 | -0.003 | 1.54e-02 | Kcnj2 | -0.541 | 3.30e-04 | Cc2d2a | -0.345 | 1.53e-04 |
| 5 | Sema6b | 0.357 | 1.78e-05 | Ppp1r12b | 0.191 | 6.12e-02 | Ndfip2 | -0.29 | 5.50e-04 | Zbtb16 | 0.81 | 2.52e-04 |
| 6 | Ogt | -0.223 | 1.78e-05 | Smc6 | -0.435 | 8.44e-02 | Sema4b | 0.383 | 6.85e-04 | Mpdz | 0.453 | 5.83e-04 |
| 7 | Cyld | -0.256 | 1.78e-05 | Comtd1 | 0.676 | 8.44e-02 | Tmem201 | 0.219 | 8.64e-04 | Hspa5 | 0.364 | 2.58e-03 |
| 8 | Hnrnpa0 | 0.327 | 6.34e-05 | Tnfrsf11b | -0.011 | 8.44e-02 | Prrt1 | 0.335 | 2.01e-03 | Cables1 | 0.443 | 2.95e-03 |
| 9 | Etnk1 | -0.165 | 3.52e-04 | Nr1d2 | -0.192 | 1.06e-01 | Grin2d | 0.483 | 2.39e-03 | Rere | 0.268 | 3.07e-03 |
| 10 | Derpc | -0.006 | 3.52e-04 | Plk4 | -0.012 | 1.06e-01 | Sulf2 | 0.284 | 2.39e-03 | Zkscan2 | -0.393 | 4.29e-03 |

### Heatmap of Top Differentially Expressed Genes
A heatmap of the top 50 most significant DEGs from the CR vs Chow comparison was generated from VST-normalized counts, with rows scaled by z-score to enable cross-gene comparison (Figure 2). Two distinct gene clusters were visible: one showing coordinated upregulation in CR samples relative to all other dietary groups, and one showing coordinated downregulation. The CR group displayed the most distinct expression pattern, consistent with its status as the comparison group for gene selection. Partial expression of the CR-associated pattern was also observable in the 5% protein group, with several genes showing intermediate expression levels relative to Chow and CR, while the 10% and 15% protein groups showed expression patterns more similar to Chow. Sex-associated variation was visible in the column annotations, consistent with the transcriptional sex differences identified in the sample distance heatmap.

<img width="1676" height="1245" alt="image" src="https://github.com/user-attachments/assets/f3a8af1e-9a3b-4ec0-8688-8c07b6cdad93" />

&nbsp;
**Figure 2. Heatmap of the top 50 most significant differentially expressed genes from the CR vs Chow comparison.** 
Genes were selected by ranking all DESeq2 results by adjusted p-value and taking the top 50. Expression values are variance-stabilizing transformed (VST) counts, scaled by z-score across rows to enable cross-gene comparison. Columns are ordered by diet group (CR, P5, P10, P15, Chow) and annotated by diet and sex. Row clustering is based on Euclidean distance with complete linkage.

### Volcano Plot — Caloric Restriction vs Chow
Differential expression results for the CR vs Chow comparison were visualised using a volcano plot displaying log2 fold change against -log10 p-value for all 15,133 tested genes (Figure 3). Of 366 significant DEGs (padj < 0.05), 194 were upregulated and 172 were downregulated in CR relative to Chow. The majority of significant genes showed modest fold changes consistent with apeglm shrinkage estimation, with most falling between -0.5 and 0.5 log2 fold change. The most significantly upregulated genes included Dbp (log2FC = 0.58, padj = 1.27×10⁻⁵), Cited2 (log2FC = 0.48, padj = 6.93×10⁻⁶), and Sema6b (log2FC = 0.36, padj = 1.78×10⁻⁵), while the most significantly downregulated genes included Bmal1 (log2FC = -0.40, padj = 5.84×10⁻⁶) and Clock (log2FC = -0.25, padj = 6.93×10⁻⁶).

<img width="2700" height="2100" alt="image" src="https://github.com/user-attachments/assets/bf389c47-de6f-43ad-8b3b-f5f2b0effef4" />

&nbsp;
**Figure 3. Volcano plot of differential gene expression in the CR vs Chow comparison.** 
Each point represents one gene, plotted by log2 fold change (x-axis) against -log10 p-value (y-axis). Significant DEGs (padj < 0.05) are highlighted in red (upregulated) or blue (downregulated). Non-significant genes are shown in grey. The top 15 most significant DEGs are labelled by gene symbol. Fold change estimates are apeglm-shrunken. Dashed vertical lines indicate log2 fold change of ±1; dashed horizontal line indicates padj = 0.05.

### Gene Set Enrichment Analysis
Gene set enrichment analysis (GSEA) was performed for all four dietary comparisons using genes ranked by DESeq2 Wald test statistic, allowing detection of coordinated pathway-level changes across the full transcriptome independent of significance thresholds. The number of significantly enriched GO Biological Process terms (nominal p-value < 0.05) varied across comparisons: CR yielded 31 enriched terms, 10% protein yielded 26, 15% protein yielded 65, and 5% protein yielded 8.

In the CR comparison, activated pathways included glycolytic and glucose metabolic processes, glial cell differentiation, and regulation of neurogenesis, while sensory perception pathways were suppressed (Figure 4a). In the 10% protein comparison, canonical glycolysis pathways were activated alongside suppression of circadian rhythm and regulation of synaptic vesicle exocytosis pathways (Figure 4b). The 15% protein comparison was dominated by activation of protein folding, endoplasmic reticulum stress response, and chaperone-mediated protein folding pathways, with no suppressed terms reaching significance (Figure 4c). The 5% protein comparison yielded 8 suppressed terms including cilium movement and potassium ion transport, with no activated terms reaching significance (Figure 4d).

<img width="1100" height="900" alt="GSEA_CR" src="https://github.com/user-attachments/assets/9994f86c-6b3f-4aba-8544-aaa7c5d75ed1" />
<img width="1100" height="900" alt="GSEA_P10" src="https://github.com/user-attachments/assets/d7cab23e-b2c1-4f8c-adaa-499f110e9036" />
<img width="1100" height="900" alt="GSEA_P15" src="https://github.com/user-attachments/assets/5add4bef-69cb-4456-b796-414bfb29b5d5" />
<img width="1100" height="900" alt="GSEA_P5" src="https://github.com/user-attachments/assets/04096c97-487e-4dec-acef-225f566093e0" />
Figure 4a. GSEA dotplot for CR vs Chow. Figure 4b. GSEA dotplot for 10% Protein vs Chow. Figure 4c. GSEA dotplot for 15% Protein vs Chow. Figure 4d. GSEA dotplot for 5% Protein vs Chow. Dotplots show the top 15 enriched GO Biological Process terms, split by activation (positive enrichment score) and suppression (negative enrichment score). Dot size represents the number of genes in the gene set; dot colour represents adjusted p-value. Genes were ranked by DESeq2 Wald test statistic. Gene sets between 15 and 500 genes were considered.

---

## Discussion 
This analysis re-examined hippocampal RNA-seq data from Wahl et al. (2018) using a modern pseudoalignment-based quantification pipeline and broadly replicated the paper's key transcriptional findings while revealing additional pathway-level insights through GSEA. Caloric restriction produced the most pronounced transcriptional response (366 DEGs), followed by 10% protein (116 DEGs), 15% protein (52 DEGs), and 5% protein (4 DEGs), consistent with the original study's finding that CR induces more extensive hippocampal gene expression changes than LPHC diets. The non-monotonic relationship between protein restriction severity and DEG count likely reflects high within-group biological variability at dietary extremes rather than a true absence of transcriptional response, a limitation of the small sample size (n=6 per group).

### Circadian Rhythm Remodelling Under Caloric Restriction
The most interesting finding of this analysis was the disruption of circadian rhythm gene expression under CR, with Bmal1 and Clock downregulated and Dbp upregulated as top CR-associated DEGs. BMAL1 and CLOCK form a heterodimeric transcription factor complex that drives expression of clock-controlled output genes including Dbp, while period and cryptochrome proteins feedback to suppress their activity (Takahashi, 2017). The opposing directional changes observed here suggest CR remodels the phase and amplitude of circadian oscillations in hippocampal tissue rather than simply suppressing or activating the clock. This is biologically significant as hippocampal circadian rhythms regulate memory consolidation, synaptic plasticity, and cognitive aging (Smarr et al., 2014), and CR-induced remodelling of these rhythms may contribute to the cognitive benefits of dietary restriction. That GSEA did not identify circadian rhythm as enriched under CR likely reflects signal concentration in a small number of highly significant individual genes rather than a distributed gene set effect, illustrating a known limitation of set-level testing methods.

### Gpr17 as a Shared Transcriptional Response to Protein Restriction
Gpr17 was the most significantly upregulated gene across all three protein restriction groups, directly replicating Wahl et al. (2018). Gpr17 encodes an orphan G-protein coupled receptor expressed in oligodendrocyte precursor cells that regulates the transition from immature to mature myelinating oligodendrocytes (Chen et al., 2009). Its consistent upregulation across all degrees of protein restriction but not under CR suggests that reduced dietary protein specifically activates oligodendrocyte precursor signalling independent of caloric intake. This is supported by the GSEA finding that glial cell differentiation was activated under both CR and P15, though the shared involvement of Gpr17 in neuroinflammatory responses means its upregulation could alternatively reflect a nutritional stress response rather than a beneficial adaptive change (Lewash et al., 2025).

### Pathway-Level Insights and Methodological Considerations
GSEA revealed coordinated pathway-level changes that were not captured by individual gene analysis alone. Under CR, the activation of glycolytic and glucose metabolic pathways is consistent with the well-established role of caloric restriction in shifting cellular metabolism toward increased glucose utilisation (Fontana & Partridge, 2015). The 10% protein comparison showed a distinct pattern from CR, with suppression of circadian rhythm and synaptic signalling pathways suggesting that moderate protein restriction and caloric restriction affect hippocampal gene expression through partially distinct mechanisms. The 15% protein comparison was dominated by activation of protein folding and ER stress pathways, consistent with Hspa5 as a top DEG which suggest that protein restriction at this level engages a cellular stress response in hippocampal tissue (Walter & Ron, 2011). 

Methodologically, despite differences between the pipeline used here and that of Wahl et al. (who used TopHat2/featureCounts with FPKM normalisation compared to the Salmon/DESeq2 approach used in this analysis) strong concordance in top DEGs across both studies suggests the core biological findings are robust. The inclusion of sex as a covariate, justified by the sex-based clustering observed in the sample distance heatmap and the sex-specific protein expression differences reported by Wahl et al., represents a more statistically rigorous approach that improved power to detect diet-specific transcriptional effects.

Taken together, these findings suggest that caloric restriction and dietary protein restriction engage partially overlapping but still distinct transcriptional mechanisms in the aging hippocampus, with Gpr17 emerging as a consistent molecular signature of protein restriction and circadian rhythm remodelling as a hallmark of caloric restriction, highlighting the potential for diet to modulate hippocampal gene expression in ways relevant to brain aging and cognitive health.


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
- Takahashi JS. (2017). Transcriptional architecture of the mammalian circadian clock. Nature Reviews Genetics, 18(3), 164–179.
- Smarr BL, et al. (2014). A time to remember: the role of circadian clocks in learning and memory. Behavioural Neuroscience, 128(3), 283–303.
- Chen Y, et al. (2009). The oligodendrocyte-specific G protein-coupled receptor GPR17 is a cell-intrinsic timer of myelination. Nature Neuroscience, 12, 1398–1406.
- Lewash M, et al. (2025). orphan G protein-coupled receptor with therapeutic potential. Trends in Pharmacological Sciences, Cell Press.
- Fontana L, Partridge L. (2015). Promoting health and longevity through diet: from model organisms to humans. Cell, 161(1), 106–118.
- Walter P, Ron D. (2011). The unfolded protein response: from stress pathway to homeostatic regulation. Science, 334(6059), 1081–1086.

