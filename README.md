## chip_seq_snakemake  

This is a simple utility to analyse CHIP-seq data with a selection of popular tools. It is written in `snakemake` that allows the user to scale it to any system. Currently, it contains analyses with the following tools:  

+ `fastqc`(quality control for reads)  
+ `trim_galore`(Adapter trimming)  
+ `cutadapt`(Adapter trimming)  
+ `bowtie2`(alignment)  
+ `picard`(Mark duplicates)  
+ `deeptools`(Generate gene-body meta-profile from bigWig files, Calculate genome-wide IP enrichment relative to control)  
+ `phantompeakqualtools`(Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC)  
+ `MACS2`(Call broad/narrow peaks)  
+ `HOMER`(Annotate peaks relative to gene features)  

### Installation  
For installation purposes, it is a good idea to make a new python virtual environment and install snakemake on it, and then download all the tool binaries in your desired location. Following this, retrieve the repository using `git`;  

```sh
git clone https://github.com/Moonerss/chip_seq_snakemake.git
```

The only other modification that needs to be done is to edit the `config/config.yaml` and `config/samplesheet_test.csv` files to specify the location of files on your system, and then run following in the repository folder  

