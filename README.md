# 🏥 Nile International Hospital - Cumulative Antibiogram

> **Source of Health**  
> *Developed by Aslam Kimbugwe - Medical Laboratory Scientist*

![Project Status](https://img.shields.io/badge/Status-Live-success)
![Python](https://img.shields.io/badge/Python-3.x-blue?logo=python)
![Plotly](https://img.shields.io/badge/Visuals-Plotly-cyan?logo=plotly)

## 📖 Overview
This project provides a "World-Class" automated system for generating **Cumulative Antibiograms** based on Clinical and Laboratory Standards Institute (**CLSI M39**) guidelines. It transforms raw laboratory AST (Antimicrobial Susceptibility Testing) data into interactive, branded, and actionable reports for antimicrobial stewardship.

## ✨ Key Features
- **CLSI M39 Compliance**: Automated deduplication logic (keeping only the first isolate per patient/species per period).
- **Interactive Dashboard**: High-quality HTML heatmaps with hover data for exact susceptibility percentages and isolate counts ($N$).
- **Specimen Stratification**: Automatic breakdown of susceptibility trends for the Top 5 most frequent sample types (e.g., Blood, Urine, Pus Swab).
- **Gram-Specific Analysis**: Dedicated sections for Gram-Negative and Gram-Positive pathogens.
- **Fungi & Yeast Monitoring**: Integrated monitoring for Candida species filtered strictly by Antifungals.
- **Professional PDF Export**: One-click generation of high-resolution, branded PDF reports for physical distribution and presentations.
- **Branded Design**: Clean, modern interface using Nile International Hospital's official brand colors.

## 📁 Project Structure
- `generate_antibiogram.py`: The core engine for data processing and report generation.
- `Antibiogram_Report.html`: The interactive web dashboard.
- `Antibiogram_Report.pdf`: The print-ready high-resolution report.
- `.gitignore`: Configured to protect sensitive patient information (PII/PHI).

## 🚀 Getting Started

### Prerequisites
- Python 3.8+
- Required Libraries:
  ```bash
  pip install pandas plotly kaleido fpdf2
  ```

### Usage
1. Place your laboratory data files in the `Data/` folder.
   - `Isolates per specimen.csv`
   - `AST.csv`
2. Run the generator script:
   ```bash
   python generate_antibiogram.py
   ```
3. Open `Antibiogram_Report.html` in your browser or share the generated `Antibiogram_Report.pdf`.

## 🛡️ Data Privacy
This repository is configured to **strictly exclude** patient-identifiable data. Only the source code and report templates are tracked. Raw data files (`.csv`) are ignored to maintain medical data confidentiality.

---
*© 2026 Nile International Hospital. Confidential Internal Use Only.*
