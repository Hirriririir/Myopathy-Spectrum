# Myopathy spectrum
> A integration analysis of bulk RNA-seq data from human skeletal muscles (1221 muscles x 9231 genes). 

The final integration dataset is shared in [h5ad](https://anndata.readthedocs.io/en/latest/) format (a format commonly used for single-cell dataset), which has the advantage of including multiple layers of count and meta data in one single file. 


## Inclusion and exclusion criteria
- Only human skeletal muscle tissue (no cell lines or organoids) 
- Bulk-RNA sequencing by high throughput technics (no chip arrays)
- Raw count data preserved (datasets shared in transformed count format were excluded)



## Folder descriptions
- **DEG**: Differential expression analysis (DEG) results exported from edgeR.
- **Meta**: All meta data of the integration dataset (Helsinki/GEO/GTEx): table columns including sample_id, phenotype, geo_accession (batch), sequencing method, platform, age range, gender.
- **Helsinki data**: 127 skeletal muscle bulk RNA-seq data from Helsinki ([Group Udd](https://www.folkhalsan.fi/en/knowledge/research/genetics/group-udd/), Folkhälsan Research Center, University of Helsinki).
- **GEO data**: 291 skeletal muscle bulk RNA-seq data downloaded from the GEO database ([GSE115650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115650), [GSE140261](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140261), [GSE175861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175861), [GSE184951](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184951), [GSE201255](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201255), [GSE202745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202745)).
- **GTEx data**: 803 skeletal muscle bulk RNA-seq data downloaded from the GTEx Analysis V8 ([dbGaP Accession phs000424.v8.p2](https://gtexportal.org/home/datasets#datasetDiv1)).
- **Validation**: validation data downloaded from the supplementary files from the used GEO datasets or generated from the integration dataset.

## Mainly used packages and tools
- **Python (3.8.1)**: [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) (high-dimensional data processing), [gseapy](https://gseapy.readthedocs.io/en/latest/gseapy_example.html) (pathway analysis), [TAPE](https://github.com/poseidonchan/TAPE) (celltype deconvolution), [conorm](https://gitlab.com/georgy.m/conorm) (count normalization).
- **R (4.2.2)**: [EdgeR](https://bioconductor.org/packages/edgeR/) (DEG analysis), [ComBat-seq](https://github.com/zhangyuqing/ComBat-seq) (batch adjustment).
