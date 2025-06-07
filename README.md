# DeepChek® NGS Library Preparation Script for MGISP-100

Automation script for performing DeepChek® NGS Library Preparation with UDI Adapters (MGI) on the MGISP-100 robotic system. This repository contains the validated protocol, software requirements, labware setup, and consumable information necessary to perform DNA library construction suitable for sequencing on MGI platforms (E25, G99, G400).

---

## 📖 Description

This project automates the full workflow of the **DeepChek® NGS Library Preparation (RUO)** using the **MGISP-100** system. The script replicates manual processes including DNA fragmentation, end repair, A-tailing, adapter ligation, and library amplification — turning up to 16 DNA samples into sequencing-ready libraries compatible with MGI’s DNBSEQ platforms.

The automation is adapted for DeepChek reagent kits **REF 203A24 / 203A96**, supporting both 24 and 96 reaction formats. It also integrates QC checkpoints, reagent volume calculators, and follows lab safety and equipment preparation guidelines aligned with MGI and ABL protocols.

---

## 🚀 Features

- 📦 Supports DeepChek® NGS Library Prep 24 & 96 RXN kits  
- 🧪 Fully automated for MGISP-100 with PCR, Magnet, and Temp modules  
- 🧬 Up to 16 samples per run  
- 📊 Includes library quantification and QC recommendations  
- ✅ Compatible with low/high complexity DNA and FFPE  
- 🔄 Fragmentation to ready-to-sequence libraries in a single run

---

## 🧰 Requirements

### 🔧 Hardware Modules

| Module                     | Position     |
|---------------------------|--------------|
| PCR Module                | Pos3         |
| Magnet Module             | Pos6         |
| Temp. Control Module      | Pos5 Col 6-8 |
| 8-Tube Strip Cover Module | Pos1         |
| Trash Can                 | Pos7         |

### 🧪 Reagent Kits

- DeepChek® NGS Library Preparation (REF 203A24 / 203A96)
- DeepChek® NGS Clean-up beads (REF N411-03 / N411-04)
- Adapter DNBSEQ (MGI)

### 🧫 Required Equipment

- MGISP-100 v1.9.3.476 or higher
- Qubit® 3.0 Fluorometer / Bioanalyzer / Fragment Analyzer
- Standard lab pipettes, vortex, centrifuge

### 📦 Customer-Prepared Consumables

| Consumable                                         | Brand | Cat. No.     | Quantity     |
|----------------------------------------------------|-------|--------------|--------------|
| 250 μL automated filter tips                       | MGI   | 1000000723   | 4 Boxes      |
| 1.3 mL 96-Well U-bottom Deepwell Plate             | MGI   | 1000004644   | 2 Plates     |
| Hard-shell Thin-wall 96-well Skirted PCR Plates    | MGI   | 091-000165-00| 2 Plates     |
| Break-away PCR Plate and Cover, 96-Well            | MGI   | 100-000016-00| 10 Strips    |
| 2 mL SC micro tube, PCR-PT                         | MGI   | 1000001553   | 5 Tubes      |
| 0.5 mL SC micro tube, PCR-PT                       | MGI   | 1000001558   | 2 Tubes      |


### 📂 Required Files

- Script files for MGISP-100
- `DeepChek NGS Reagent and Consumables Calculation.xlsx` (for setup)

---

## ⚠️ Precautions

- Ensure PCR program is correctly installed on the MGISP-100 before running
- Use only the recommended reagents and consumables
- Perform **Pre-clean** and **Post-clean** procedures
- Refer to the full SOP for QC and normalization guidance
- Dispose of samples and waste responsibly
- For technical support, contact: [MGI-service@mgi-tech.com](mailto:MGI-service@mgi-tech.com)

---

## 🔬 Workflow Summary

1. DNA Input (50 ng in 35 µL of molecular-grade water)
2. Fragmentation, End Repair & A-tailing
3. Adapter Ligation
4. Amplification
5. Bead Cleanup
6. QC (Optional: Bioanalyzer or Qubit)
7. Library Normalization
8. Ready for DNB Circularization and Sequencing

---

## 🌐 Website

[MGI-tech.eu](https://www.mgi-tech.eu)

---

## 🏷️ Topics

- MGISP100 Scripts
- NGS Automation
- DeepChek Library Prep
- Liquid Handling
- DNA Sequencing Preparation

---

## 📁 Repository Contents
📁 en-eu/ # Main MGISP-100 automation scripts
📁 imgs_workflow/ # Visual setup and workflow diagrams
📄 README.md # This documentation
📄 DeepChek_Calculation.xlsx # Reagent and consumable calculator
📄 Protocol_manuals.pdf # Step-by-step workflow summary
