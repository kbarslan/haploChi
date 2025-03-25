# HaploChi 🧬

> **HaploChi** – Haplotype Chimera Inspector  

---

## 🧬 What is HaploChi?

**HaploChi** is a modular pipeline for assessing haplotype resolution quality, identifying chimeric contigs, and summarizing genome assembly outputs using long and short reads.

It automates mapping, chimera detection, GFF annotation refinement, and high-level reporting – all in a single flexible workflow.

---

## 🖼️ Pipeline Overview

<p align="center">
  <img src="docs/haplochi_pipeline.png" alt="HaploChi Pipeline Overview" width="700"/>
</p>

---

## 🧠 Algorithm Explanation

<p align="center">
  <img src="docs/haplochi_pipeline_alg.png" alt="Algorithm Explanation" width="700"/>
</p>

---

## 🧪 Example Outputs

<p align="center">
  <img src="docs/vertigo.png" alt="Vertigo Example" width="340"/>
  <img src="docs/sorg_49.png" alt="Sorg_49 Example" width="340"/>
</p>

---

## 🚀 Quick Start

```bash
git clone https://github.com/yourusername/haplochi.git
cd haplochi

# Install dependencies (Python + tools)
pip install -r requirements.txt

# Run example steps
python haplochi.py split ds.list
python haplochi.py LRmapping sorg_49_HR.fa ds.list
python haplochi.py SRmapping sorg_49_HR.fa s.list
