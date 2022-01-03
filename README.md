# SRAstats
### overview statistics of the NCBI SRA database


*How to use*:
1. download RunInfo table from NCBI SRA (for example, https://www.ncbi.nlm.nih.gov/sra/?term=insecta, and then go to Send to > file > RunInfo)
2. Run SRAstats, for example:

```
Rscript --no-restore SRAstats.R -i insecta.csv -o results.dir..insecta -n Insecta
```
3. Find the generated plots and project-centric table in the output dir
