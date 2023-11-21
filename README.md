# BIOINFO_PROJECTS

🧬 **FastQC Report for SRR26670608_1.fastq**

[![Linux code](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/LINUX%20CODE.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/LINUX%20CODE.PNG)

## Overview

🔍 Welcome to the in-depth FastQC report for SRR26670608_1.fastq. This analysis was conducted using FastQC version 0.11.9, a powerful tool designed to evaluate the quality of high-throughput sequencing data.

[![Summary](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/fastqc%20Report.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/summary.PNG)

## Summary

- **Basic Statistics:** PASS [![Statistic](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/statistics.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/statistics.PNG)
  - The file SRR26670608_1.fastq is encoded with Sanger/Illumina 1.9 standards, comprising 15,884,288 sequences, each 150 base pairs long, with a balanced GC content of 50%.

- **Per Base Sequence Quality:** PASS [![Per base quality graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/2.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/2.PNG)
  - The per base quality graph illustrates consistent and high-quality sequencing across all positions.

- **Per Tile Sequence Quality:** PASS [![Per tile quality graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20tile%20sequence%20quality.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20tile%20sequence%20quality.PNG)
  - The per tile quality graph demonstrates a uniform distribution of quality scores across different tiles, ensuring reliable data acquisition.

- **Per Sequence Quality Scores:** PASS [![Per Sequence quality graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20sequence%20quality%20scores.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20sequence%20quality%20scores.PNG)
  - The per sequence quality scores graph indicates consistently high-quality scores for the entire dataset.

- **Per Base Sequence Content:** FAIL [![Per base sequence content graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20base%20sequence%20content.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20base%20sequence%20content.PNG)
  - The per base sequence content graph reveals variations, suggesting potential issues in certain regions. Further investigation is recommended.

- **Per Sequence GC Content:** PASS [![Per sequence GC content graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20sequences%20GC%20content.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20sequences%20GC%20content.PNG)
  - The per sequence GC content graph shows a balanced distribution of GC content across all sequences, indicating robust sequencing performance.

- **Per Base N Content:** PASS [![N content graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20base%20N.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/per%20base%20N.PNG)
  - The N content graph indicates a minimal presence of ambiguous bases, reflecting the overall high quality of the sequencing data.

- **Sequence Length Distribution:** PASS [![Sequence length distribution graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/sequence%20Length%20Distribution.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/sequence%20Length%20Distribution.PNG)
  - The sequence length distribution graph showcases a consistent length of 150 base pairs across all sequences.

- **Sequence Duplication Levels:** FAIL [![Duplication level graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/sequence%20duplication%20levels.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/sequence%20duplication%20levels.PNG)
  - The duplication level graph indicates potential concerns regarding sequence duplication. Further exploration is recommended.

- **Overrepresented Sequences:** PASS 
  - No overrepresented sequences were detected in the dataset.

- **Adapter Content:** PASS [![Adapter content graph](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/adapter%20content.PNG)](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/fastq_report/adapter%20content.PNG)
  - The adapter content graph confirms the absence of adapter contamination, reaffirming the data's purity and suitability for downstream analyses.

## Conclusion

The overall quality of the sequencing data in SRR26670608_1.fastq is satisfactory, with most metrics passing the quality checks. However, attention is needed for the Per base sequence content and Sequence Duplication Levels, where the analysis indicates potential issues that may require further investigation.

---

### About FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular bioinformatics tool developed by the Babraham Institute. It is designed to provide a quick overview of the quality of high-throughput sequencing data, offering valuable insights for researchers and bioinformaticians in the preprocessing of sequencing datasets.

🚀 **Explore More:** If you want to dive deeper into the results or explore other bioinformatics projects, check out the [full repository](https://github.com/Rachel2705/BIOINFO_PROJECTS). Feel free to connect and collaborate! 🧑‍💻

---

# Bioinformatics Data Visualization

Absolutely! Here's an enhanced version of your README with added icons and more visual appeal:

```markdown
# 🌐 Bioinformatics Data Visualization

Welcome to the Bioinformatics Data Visualization repository! This project involves creating an intriguing data visualization of combined sine and cosine waves using Matplotlib. The generated plot showcases three subplots, each depicting different aspects of the data.

## 📋 Requirements

- Python (>=3.6)
- ![Matplotlib](https://img.shields.io/badge/Matplotlib-v3.4.2-blue)
- ![NumPy](https://img.shields.io/badge/NumPy-v1.21.0-green)

Install the required dependencies using the following command:

```bash
pip install -r requirements.txt
```

## 🚀 Usage

Explore the data visualization script in the [data_visualization](https://github.com/Rachel2705/BIOINFO_PROJECTS/tree/main/data_visualization) directory of the [Rachel2705/BIOINFO_PROJECTS](https://github.com/Rachel2705/BIOINFO_PROJECTS) repository.

## 📊 Plot Description

The generated plot consists of three subplots:

1. **Sine Wave**: Displaying a sinusoidal function in 🟦 blue.
2. **Cosine Wave**: Displaying a cosinusoidal function in 🟩 green.
3. **Sum of Sine and Cosine**: Displaying both the sine and cosine waves along with their sum in 🟥 dashed blue, 🟥 dashed green, and 🟥 solid red lines, respectively.

The entire figure is titled "Multiple Plots of Sine and Cosine Waves."

## 🌈 Additional Features

## 🖼️ Preview


![Multiple Plots](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/data_visualization/multiple_plot.png)

![single plot](https://github.com/Rachel2705/BIOINFO_PROJECTS/blob/main/data_visualization/plot.png)

## 🛠️ Customization

Feel free to customize the script to suit your specific bioinformatics data and visualization needs.

## 🔗 Related Projects

Check out [Rachel2705/BIOINFO_PROJECTS](https://github.com/Rachel2705/BIOINFO_PROJECTS) for more bioinformatics projects, including data visualization.

## 📄 License

This project is licensed under the [MIT License](LICENSE).
```

This version includes badges for Matplotlib and NumPy versions, emojis for visual appeal, and shields for dependencies. Feel free to further customize it based on your preferences!
