# Processing input data

This project contains a small snakemake pipeline for taking BAM files to SAF files appropriate for analysing in the `window-results` project. It is used for both the human and the impala data in the manuscript.

## Getting the BAM files

### Human 1000g data

For the human data, the BAMs are available through the 1000G projects [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/). We used IDs `NA12341`, `NA07357`, `NA11933`, `NA12760`, `NA12872`, `NA12546`, `NA12827`, `NA12275`, `NA12829`, `NA07056` for CEU and `NA18877`, `NA19095`, `NA19160`, `NA18923`, `NA19130`, `NA19114`, `NA19200`, `NA19138`, `NA18864`, `NA18915` for YRI.

## Impala

For the impala data, the FASTQs are available via the SRA with accession `PRJNA862915`. Sample names `1508`, `1509`, `1510`, `1511`, `1514`, `1515`, `1516`, and `1517` are included in the Shangani population in the manuscript, while sample names `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`, `11`, `8716`, and `4009` belong to Maasai Mara. The FASTQs were cleaned and aligned by Genis, who also prepared the filters used here. Code for these steps can be found [TODO: add link].

## Running

Once BAM files are obtained, the snakemake pipeline runs from a config file each for human and impala data. Examples are included in the `configs` subdir and assume that a reference is available together with a list of chromosomes to use. In the manuscript, we used all autosomal chromosomes greater than 100kb. For the impala in particular, a sites filter was moreover used: see above for directions to how this was made.

After filling out the config files for both impala and human data, a `run.fish` script is available for shelling out all appropriate commands to snakemake. 

```shell
fish run.fish $threads
```

If you do not have fish, copying out the commands to your preferred shell should be straightforward.

By default, results will appear in the `results/` directory, where they can be used for input to the `window-results` workflow.