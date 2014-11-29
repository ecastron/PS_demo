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

### Installing PathoScope and PS_demo
PathoScope is now hosted in GitHub so you can easily get it by issuing the following command from the Terminal  

		git clone https://github.com/PathoScope/PathoScope.git

And PS_demo

		https://github.com/ecastron/PS_demo.git

### Get data and get reference genomes
We are going to use data from a study exploring microbiome diversity in oropharingeal swabs from schizophrenia patients and healthy controls. The SRA accession number is `SRR1519057`. 

![SRA](https://github.com/ecastron/PS_demo/raw/master/img01.png)

This file is probably too big for a demo so I randomly subsampled the reads down to a more manageable size (~40 M to 40 K reads)  
- Go ahead and download the data [here](https://raw.githubusercontent.com/ecastron/PS_demo/master/ES_211.fastq)
- Now you need at least two files, one as target library (where your reads are going to be mapped) and as filter library (internal controls, host genome, contaminants, etc. that you want to remove)

As target library, you can use any multi fasta file containing full or draft genomes, or even nucleotide entries from NCBI. The only condition is that the fasta entries start with the taxonomy ID from NCBI as follows:

From:  \>gi|40555938|ref|NC_005309.1| Canarypox virus, complete genome  To:  \>**ti|44088|**gi|40555938|ref|NC_005309.1| Canarypox virus, complete genome  

You could do this very easily in **PathoLib**:

		python pathoscope.py LIB -genomeFile my_file.fasta -outPrefix target_library

Alternatively, we provide the entire NCBI nucleotide database already formatted [here] (ftp://pathoscope.bumc.bu.edu/data/nt_ti.fa.gz) (10 GB file). You could also use **PathoLib** to subsample this big file (50 GB uncompressed) and select only the taxa that you want. For instance, obtaining all the virus entries in nt_ti.fa (virus taxonomy ID = 10239)

		python pathoscope.py -LIB python pathoscope.py LIB -genomeFile nt_ti.fa -taxonIds 10239 --subTax -outPrefix virus

Or in order to create a filter library, say all human sequences:
		
		python  pathoscope.py -LIB python pathoscope.py LIB -genomeFile nt_ti.fa -taxonIds 9606 --subTax -outPrefix human

However, I'm providing a target and filter library already formatted that you can download [here](https://www.dropbox.com/sh/3kjvec5lizwmo9l/AAAQmHEAwAfDixGtKC6eTeN1a?dl=0) and [here](). The target library is a collection of genomes from the reference library of the Human Microbiome Project (description [here](http://hmpdacc.org/HMREFG/)), and the filter library is simply the human genome (hg19). We are also going to use another filter library (phix174) to get rid of all the reads mapping to the Illumina internal control sequence that is sometimes added to sequencing experiments.

### Let's map the reads
Once you have your data and target and filter libraries, we are ready to go ahead with the mapping step. For this we use bowtie2 so we will need to tell **PathoMap** where the bowtie2 indices are. If you don't have bowtie2 indices, not a problem, **PathoMap** will create them for you. And if your fasta files are larger than 4.6 GB, not a problem either, **PathoMap** will split your fasta files and create indices for each one of the resulting files.

If you have fasta files and not bowtie2 indices:

		python pathoscope.py -MAP -U ES_211.fastq -targetRefFiles HMP_ref_ti_0.fa,HMP_ref_ti_1.fa -filterRefFiles human.fa,phix174.fa  -outDir . -outAlign ES_211.sam  -expTag DAV_demo

But if you already have Bowtie2 indices (our case), you can issue the following command:

		python pathoscope.py -MAP -U ES_211.fastq -targetIndexPrefixes HMP_ref_ti_0,HMP_ref_ti_1 -filterIndexPrefixes genome,phix174  -outDir . -outAlign ES_211.sam  -expTag DAV_demo

Let's give it a try...



