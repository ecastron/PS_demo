# Analyzing Read Counts from PathoScope

##### By now you already obtained read counts from your metagenomic experiment. The thing is what do we do with them. Most of the time, obtaining read counts is only the beggining to exploring the metagenomic content of your samples. The first step in analyzing read count data is to put it in a contingency table format.

### Let's combine multiple PathoScope tsv files into a contingency table

We need to create a table where rows are genomes and columns are samples. In the previous demo we analyzed a downsampled version of sample ES_211. The original dataset consisted of 16 cases (patients with schisophrenia) and 16 healthy controls.  
I've already analyzed these data and I'm providing all tsv files. If you want to give it a try, you can always download the original data which was deposited in NCBI under BioProject [PRJNA255439] (http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA255439).  
The following code is from an R script that I use but that it's not as tidy and neat as I'd like to.

Let's load the required libraries and set our working directory

```{r}
library(xlsx)
library(gtools)

directory<-"~/PS_demo/tsv"
setwd(directory)
```

If you don't have those libraries installed, simple issue the following commands

```{r}
install.packages(c("xlsx","gtools"))
```

We need to create a list of files, a list of "patient names", and finally read in relevant columns in the tsv's so that we can combine everything into a contingency table  

```{r}
# make a list of the files to loop through

list_of_files <- list.files(path=directory, pattern = ".tsv")

# make a list of just the patient names

patient_names<-NULL
for( i in 1:length(list_of_files)){
  patient_names[i]<-substring(list_of_files[i], 1, 6)
}

# read in each table

read_counts <- lapply(list_of_files, read.table, sep="\t", header = FALSE, skip =2)
read_counts <- lapply(read_counts, function(x) x[, c(1,4)])
read_counts <- lapply(read_counts, function(x) x[complete.cases(x),])
```

Now we have a list of tables (or data frames) that we need to combine, in a non-redundant fashion, to create our contingency table

```{r}
# for each table make the first col name OTU and the second the patient name

for( i in 1:length(list_of_files)){
  colnames(read_counts[[i]])<- c("OTU", patient_names[i])
}

# list of lists called otu which stores the first column otu names for each dataframe

otu<-NULL
for( i in 1:length(list_of_files)){
  otu[i]<- list(as.character(read_counts[[i]][, 1]))
}

# for each dataframe in read_counts transpose and then 

read_counts <- lapply(read_counts, function(x) t(x[,2]))

# add the otus back as the column name

for( i in 1:length(list_of_files)){
  read_counts[[i]]<-data.frame(read_counts[[i]])
  colnames(read_counts[[i]])<-otu[[i]]
  read_counts[[i]]<-data.frame(patient = patient_names[i], read_counts[[i]])
}

# combine the different dataframes together

otu_table <- read_counts[[1]]
for( i in 2:length(list_of_files)){
  otu_table <- smartbind(otu_table, read_counts[[i]], fill = 0)
}

# transpose the table back so that the microbes are the rows and the patients are the col

otu_table<-t(data.matrix(otu_table))
colnames(otu_table)<-patient_names
otu_table<-otu_table[2:nrow(otu_table), ]

# remove zeroes
otu_table_noZeroes<-otu_table[apply(otu_table[,-1], 1, function(x) !all(x==0)),]

write.csv(otu_table_noZeroes,"otu_table.csv")
```