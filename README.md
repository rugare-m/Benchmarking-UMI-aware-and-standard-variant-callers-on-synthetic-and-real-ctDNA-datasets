# Benchmarking UMI-aware and Standard Variant Callers on Synthetic and Real ctDNA Datasets

### Background
This project was part of my PhD program at the University of Liverpool's Health Data Science between 2020 and 2024. A preprint is currently under review with BMC Genomics and [available on Research Square.](https://www.researchsquare.com/article/rs-3610989/v1) 

### Preprint Abstract
Circulating tumour DNA (ctDNA) is a subset of cell free DNA (cfDNA) released by tumour cells into the bloodstream. Circulating tumour DNA has shown great potential as a biomarker to inform treatment in cancer patients. Collecting ctDNA is minimally invasive and reflects the entire genetic makeup of a patient’s cancer. ctDNA variants in NGS data can be difficult to distinguish from sequencing and PCR artefacts due to low abundance, particularly in the early stages of cancer. Unique Molecular Identifiers (UMIs) are short sequences ligated to the sequencing library before amplification. These sequences are useful for filtering out low frequency artefacts. The utility of ctDNA as a cancer biomarker depends on accurate detection of cancer variants. In this study, we benchmarked six variant calling tools, including two UMI-aware callers for their ability to call ctDNA variants. The standard variant callers tested included Mutect2, bcftools, LoFreq and FreeBayes. The UMI-aware variant callers benchmarked were UMI-VarCal and UMIErrorCorrect. We used both real and synthetic datasets, with and without UMI sequences. Variant callers displayed different preferences for sensitivity and specificity. Mutect2 showed high sensitivity, while returning more privately called variants than any other caller in non-UMI data – an indicator of false positive variant discovery. In UMI encoded data, UMI-VarCal detected fewer putative false positive variants than all other callers in synthetic datasets. UMI-VarCal also called the highest percentage of COSMIC variants in real samples, and only 4.4% uniquely called variants indicating high sensitivity and specificity. Our results indicate UMI-aware variant callers have potential to improve sensitivity and specificity in calling ctDNA variants over standard variant calling tools. There is a growing need for further development of UMI-aware variant calling tools if effective early detection methods for cancer using ctDNA samples are to be realised.

### Variant Filtering Strategy
![Chapter One Paper](https://github.com/rugare-m/Benchmarking-UMI-aware-and-standard-variant-callers-on-synthetic-and-real-ctDNA-datasets/assets/88198662/ad3314ae-117f-4109-a418-e74f08908302)
