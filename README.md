## Project Overview

This repository contains the data, analysis scripts, and results for a coordinate-based meta-analysis (CBMA) focusing on self-referential processing


## GENERAL INFORMATION

**Title of Project**  
**Shared brain basis for altered self-referential processing across psychiatric disorders?  
A systematic review and meta-analysis of neuroimaging studies**


## Computational Environment

Analyses were conducted using:

- Python 3.12.9
- R 4.5.2

## Repository Structure

```text
.
├── Data
│   ├── AnalysisData        # Processed data used for meta-analyses
│   ├── InputData           # Raw extracted coordinates and study-level information
│   └── neurosynth          # Neurosynth-related decoding inputs and outputs
│
├── Output
│   ├── 1_ALE               # ALE meta-analysis results
│   ├── 2_Contrast          # Contrast and comparison analyses
│   ├── 3_Tables            # Summary tables and statistics
│   ├── 4_Decoding          # Functional decoding results
│   ├── 5_Sankey            # Sankey diagrams linking disorders and brain regions
│   ├── 8_wordcloud         # Word cloud visualizations            
│   └── 9_other             # Additional or auxiliary outputs
│
├── Scripts
│   ├── 1_ALE.ipynb
│   ├── 2_Contrast.ipynb
│   ├── 3_Tables.ipynb
│   ├── 4_Decoding.ipynb
│   ├── 5_Sankey.ipynb
│   └── 6_Data_distribution.R
│
├── voxel2surface
│   ├── Conj_L.shape.gii
│   ├── ...
│   └── sz_all.png
│
├── README.md
└── requirements.txt

```

**about the research data more details can see [scientific data bank](https://www.scidb.cn/en/detail?dataSetId=1b53c91112024d71b95c42bbd748141f&version=V4).**

## Contact

### First Author

- Name: Shanshan Zhu
- ORCID: https://orcid.org/0009-0007-6612-6190
- Institution: NanJing Normal Univeristy
- Email: zhushanshan0717@gmail.com

### Corresponding Author

- Name: Hu Chuan-Peng
- ORCID: xxxxx
- Institution: NanJing Normal Univeristy
- Email: hcp4715@hotmail.com

## Contributors

- [2021.11-present] Hu Chuan-Peng: Project design; data verification and management; manuscript writing and revision.
- [2023.01-present] Shan-shan Zhu: Data collection and proofreading; data summarization and organization; manuscript writing and revision.
- [2025.03-present] Xue-Yang Zhu: Data collection and proofreading.
- [2025.03-present] Xin-Yan Li: Data collection and proofreading.
- [2025.03-present] Zhao-Li Fan: Data collection and proofreading.
- [2024.07-2024.11] Si-Yu Wu: Data collection and proofreading.
- [2023.01-2023.12] Jia-Qi Wu: Data collection and proofreading, summarization and collation.
- [2023.07-2023.09] Ya-Qi Li: Database collection and proofreading.
- [2021.11-2022.12] Shu-Ting Sun: Data collection and proofreading, summarization, data analysis, and paper writing and revision.
- [2021.11-2022.12] Nan Wang: Data collection and proofreading, summarization and collation, and paper writing and revision.
- [2021.11-2022.07] Jia-Hui Wen: Data collection.
- [2023.07-2023.11] Jian Xiao: Data collection and proofreading, database verification, data analysis.