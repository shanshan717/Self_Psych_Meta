# Shared brain basis for altered self-referential processing across psychiatric disorders?

*A systematic review and meta-analysis of neuroimaging studies*

## 📖 Project Overview

This repository contains the data, analysis scripts, and outputs for our systematic review and coordinate-based meta-analysis (CBMA) of altered self-referential processing across psychiatric disorders.

The primary meta-analysis was conducted using the Activation Likelihood Estimation (ALE) framework implemented in [NiMARE](https://nimare.readthedocs.io/).


## 🖥 Computational Environment

Analyses were conducted using:
- **Python**: 3.12.9
- **R**: 4.5.2  
- **Operating System**: macOS Sonoma 14.5

For reproducibility, Python dependencies are listed in [`requirements.txt`](requirements.txt), and R package dependencies are summarized in [`r-packages.txt`](r-packages.txt).
Because this repository does not yet include a fully locked project environment, these files should be treated as dependency manifests rather than complete lockfiles.

**Environment Setup**

```bash
# Install Python dependencies
pip install -r requirements.txt
```

Before running the R scripts in [`2_Scripts`](2_Scripts), please install the packages listed in [`r-packages.txt`](r-packages.txt).

## ✍️ Contact

If you have questions about the data or analysis pipeline, please contact:

*  **First Author**: Shanshan Zhu ([zhushanshan0717@gmail.com](mailto:zhushanshan0717@gmail.com)) | [ORCID](https://orcid.org/0009-0007-6612-6190)
*  **Corresponding Author**: Hu Chuan-Peng ([hcp4715@hotmail.com](mailto:hcp4715@hotmail.com)) | [ORCID](https://orcid.org/0000-0002-7503-5131)
*  **Institution**: Nanjing Normal University

## 📂 Repository Structure

The project is structured to ensure fully reproducible data analyses:

```text
.
├── 1_Data/                           # Input datasets for meta-analysis
│   ├── RawData/                      # Original literature data/coordinates
│   └── AnalysisData/                 # Preprocessed Sleuth format (.txt) files
├── 2_Scripts/                        # Analysis scripts (Python & R)
│   ├── 1_ALE.ipynb                   # Main ALE meta-analysis and thresholding
│   ├── 2_Contrast.ipynb              # Contrast analyses between groups
│   ├── 3_Tables.ipynb                # Generate statistical summary tables
│   ├── 4_Decoding.ipynb              # Meta-analytic functional decoding
│   ├── 5_Supply_Figure_S1_and_S2.R   # R script for supplementary figures
│   ├── 6_Diagnostics_ALE.py          # FocusCounter diagnostics
│   ├── 7_Supply_AC1.R                # Inter-rater reliability (AC1)
│   ├── 8_Supply_FSN.ipynb            # Fail-safe N (FSN) robustness analysis
│   └── 9_Supply_Visualization.ipynb  # Brain visualization
├── 3_Output/                         # Analysis results
│   ├── 1_ALE/                        # ALE z-maps and unthresholded maps
│   ├── 2_Contrast/                   # Spatial contrast maps
│   ├── 3_Tables/                     # Statistical summary tables (.tsv)
│   ├── 4_Decoding/                   # Functional decoding results & wordclouds
│   ├── 5_Supply_info_figure/         # Supplementary information figures
│   ├── 6_FocusCounter/               # Diagnostic cluster contribution data
│   ├── 7_AC1/                        # Inter-rater reliability results
│   ├── 8_FSN/                        # Fail-safe N maps and statistics
│   └── Visualization_by_Workbench/   # Neuroimaging visualization data
├── requirements.txt                  # Python dependencies
├── r-packages.txt                    # R dependencies
└── README.md
```

## 📊 Data Availability

The complete neuroimaging coordinate dataset, sample information, and unthresholded statistical maps are publicly accessible via **[Scientific Data Bank (SciDB)](https://www.scidb.cn/en/detail?dataSetId=1b53c91112024d71b95c42bbd748141f&version=V4)**.

## 👥 Contributors

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
