# PathoScope Demo for DAV
-------------------------

In this demo we will explore how to get a taxonomic profile from a metagenomic experiment using PathoScope 2.0. a more in-depth tutorial can be found in the **PathoScope** repo [here](https://github.com/PathoScope/PathoScope/raw/master/pathoscope2.0_v0.02_tutorial.pdf)

### What PathoScope can do for you
**PathoScope** is a modular piece of software that will allow you to go all the way from a fastq file to a text file (typically tab-delimited) with columns representing genomes, their proportions, etc.  
There are 6 **PathoScope modules**, however, for this demo we will focus on the three most important ones:
- ***PathoLib*** - Allows user to automatically generate custom reference genome libraries for specific scenarios or datasets
- ***PathoMap*** - Aligns reads to target reference genome library and removes sequences that align to the filter and host libraries
- ***PathoID*** - Reassigns ambiguous reads, identifies microbial strains present in the sample, and estimates proportions of reads from each genome  

Once you run your samples through **PathoScope**, you can easily import the outputfiles into R for downstream exploratory data analysis and statistical inferences.

### PathoScope Dependencies
The only dependencies for **PathoScope** are [*Bowtie2*](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [python](https://www.python.org) *2.7.3* or higher. Make sure that both are in your PATH by issuing something like `echo $PATH`

### Installing PathoScope
PathoScope is now hosted in GitHub so you can easily get it by issuing the following command from the Terminal  

		git clone https://github.com/PathoScope/PathoScope.git

### Get data and get reference genomes
We are going to use data from a study exploring microbiome diversity in oropharingeal swabs from schizophrenia patiens and healthy controls. The SRA accession number is `SRR1519057`. 

![SRA](/Users/Crandalllab/Desktop/Screen Shot 2014-11-29 at 12.37.40 PM.png)

This file is probably too big for a demo so I randomly subsampled the reads down to a more manageable size (~40 M to 40 K reads).

- Go ahead and download data here 