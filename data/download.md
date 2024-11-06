# Data
The following transcriptome annotations were used in experiments. They need to be downloaded to `data` dir and named to the corresponding name, e.g. `data/dm6.gtf`.


| Species    | Annotation               | 
| ---------- | ------------------------ | 
| Human      | GENCODE-GRCh38           | 
| Human      | RefSeq-GRCh38            | 
| Mouse      | GRCm39                   | 
| Drosophila | dm6                      | 
| Zebrafish  | GRCz11                   | 
| Maize      | Zm-B73-REFERENCE-NAM-5.0 | 
| Arabdopsis | TAIR10                   | 

# Additional assembly

To evaluate real data support for predictions by TENNIS, we downloaded a gtf assembly from a previous publication [Carlos et. al. (2023) Cell 186(11):2438-2455.e22](https://doi.org/10.1016/j.cell.2023.04.012)

```sh
mkdir GSE203583
cd GSE203583
# download file GSE203583_CIA.assembly.allTissues59K.gtf.gz
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE203583&format=file&file=GSE203583%5FCIA%2Eassembly%2EallTissues59K%2Egtf%2Egz" -O GSE203583_CIA.assembly.allTissues59K.gtf.gz
cd ..
```


