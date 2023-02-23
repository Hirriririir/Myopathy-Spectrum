# Myopathy spectrum
A integration analysis of bulk RNA-seq data from human skeletal muscles.


## Inclusion and exclusion criteria
- Only human skeletal muscle tissue (no cell lines or organoids) 
- Bulk-RNA sequencing by high throughput technics (no chip arrays)
- Raw count data preserved (datasets shared in transformed count format were excluded)




## Folder descriptions
- **DEG**: Differential expression analysis (DEG) results exported from edgeR.
- **Meta**: All meta data of the integration dataset (Helsinki/GEO/GTEx): Sample_id, phenotype, Geo_accession (batch), sequencing method, platform, age range, gender.
- **Helsinki data**: 127 skeletal muscle bulk RNA-seq data from Helsinki ([Group Udd](https://www.folkhalsan.fi/en/knowledge/research/genetics/group-udd/), Folkhälsan Research Center, University of Helsinki).
- **GEO data**: 291 skeletal muscle bulk RNA-seq data from the GEO database ([GSE115650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115650), [GSE140261](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140261), [GSE175861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175861), [GSE184951](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184951), [GSE201255](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201255), [GSE202745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202745)).
- **GTEx data**: 803 skeletal muscle bulk RNA-seq data from the GTEx Analysis V8 ([dbGaP Accession phs000424.v8.p2](https://gtexportal.org/home/datasets#datasetDiv1)).

