# Scarfo_HEC2023

This repository includes essential scripts for producing final tables and figures embedded in the following manuscript:  
**Scarfò, Randolph et al.,** _CD32 captures committed haemogenic endothelial cells during human embryonic development_ 

PubMed: 38594587\
DOI:  10.1038/s41556-024-01403-0\
GEO:  GSE223223

---

The repository is divided in three folders, 2 bulk RNAseq analyses and 1 scRNAseq analysis.  
Below a brief description of the bulk and scRNAseq workflows adopted in this work.

**Bulk RNAseq** analysis was performed using a standard pipeline that includes the follwing steps:
1. Quality control by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Trimming of bad quality reads with [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)<details><summary>Running command</summary>trim_galore --quality 20 --fastqc --length 25 --output_dir {outdir} --paired {input.r1} {inout.r2}</details>
3. Alignment with [STAR](https://github.com/alexdobin/STAR)
    <details><summary>Running command</summary>
            "STAR " +
            "--runThreadN {threads} " +
            "--genomeDir {input.genome} " +
            "--readFilesIn {params.trim_seq} " +
            "--outSAMstrandField intronMotif " +
            "--outFileNamePrefix {params.aln_seq_prefix} " +
            "--outSAMtype BAM SortedByCoordinate " +
            "--outSAMmultNmax 1 " +
            "--outFilterMismatchNmax 10 " +
            "--outReadsUnmapped Fastx " +
            "--readFilesCommand zcat "
    </details>
4. Gene expression quantification with [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889)
    <details><summary>Running command</summary>
            "featureCounts " +
            "-a {input.annot} " +
            "-o {output.fcount} " +
            "-g gene_name " +
            "-p -B -C " +
            "-s {params.strand} " +
            "--minOverlap 10 " +
            "-T {threads} " +
            "{input.bams} "
    </details>
5. Differential Expression analysis with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).  
    For Differential Gene Expression analysis we followed the standard workflow provided by package.  
   
6. Dowstream functional Analysis with [ClusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).  
    In order to retrieve functional annotation from DE analysis, we performed **O**ver **R**epresentation **A**nalysis by using the _enrichr_ function provided by the package.  
    **ORA** analysis was performed in particular using the C5 gene set from the MSigDB database (version 7.2)

---
    
**scRNAseq** analysis was performed using a standard pipeline that includes the following steps:

scRNAseq analysis was performed with [Seurat](https://satijalab.org/seurat/). Below are the main steps of the basic data analysis workflow that start from a minimal object after loading of 10X data to markers identification:  

1. Quality control and filtering
2. Cell cycle scoring
3. Normalization (default seurat settings)
4. Scaling (with following variables to regress out: percent.mt + nCount_RNA and CC.Difference calculated as show in [vignette](https://satijalab.org/seurat/articles/cell_cycle_vignette.html#alternate-workflow-1))
4. Dimensionality reduction: PCA
5. Clustering
6. Markers identification

- [6.1] Clusters related markers
- [6.2] Intracluster differential expression analysis according to comparison of interest

Pseudotime analysis was performed using [Monocle3](http://cole-trapnell-lab.github.io/monocle3/) and PAGA-tree from [dynverse](https://dynverse.org/)

In-silico perturbation analysis was performed using [CellOracle](https://github.com/morris-lab/CellOracle)

Input files for scRNAseq analysis are available in the following [link](https://www.dropbox.com/sh/83dxrxqer8cl081/AAC8zALRuYRGh1mEe4lZuaHZa?dl=0)


---

An interactive querying and exploration of the scRNAseq dataset is available at:
http://bioinfotiget.it/he/

---




