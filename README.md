# TENNIS-test

This is the test repository for [TENNIS](https://github.com/Shao-Group/TENNIS), a tool to identify missing transcripts from transcriptome annotations.

Scripts and documentation for test purposes are available in the `analyses` directory.


# Installation
Software programs should be installed and softlinked in `programs` dir. For example, TENNIS should be installed and softlinked to `programs/tennis.py`.
[cvbio](https://github.com/clintval/cvbio#cvbio) is needed to convert chromosome names between different annotations. It should be softlinked to `programs/cvbio`.  

More detailed descriptions can be found in [programs/install.md](programs/install.md)

# Data
The following transcriptome annotations were used in experiments. They need to be downloaded to `data` dir and named to the corresponding name, e.g. `data/dm6.gtf`. More about [data](data/download.md).


| Species    | Annotation               | 
| ---------- | ------------------------ | 
| Human      | GENCODE-GRCh38           | 
| Human      | RefSeq-GRCh38            | 
| Mouse      | GRCm39                   | 
| Drosophila | dm6                      | 
| Zebrafish  | GRCz11                   | 
| Maize      | Zm-B73-REFERENCE-NAM-5.0 | 
| Arabdopsis | TAIR10                   | 

# Experiments
Scripts and documentation for test purposes are available in the `analyses` directory.

Section 3.1 is reproducible using `analyze_annotations.ipynb`

Section 3.2 is reproducible using `long_read_support.ipynb`

Section 3.3 is reproducible using `pred_removed_isoforms.ipynb`