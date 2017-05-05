# SpliceVariantAnalyzer

#### Table of Contents  
[I. Introduction](#introduction)  
[II. Installation Instructions](#installation)  
[III. Usage Instructions and Examples](#instruct)   
[IV. License](#license)   
[V. Funding and Acknowledgements](#funding)   
[IV. Contact](#contact)   


<a name="introduction"/>   

## I. Introduction

SpliceVariantAnalyzer (SVA) is a splice site variant diagnosis tool which outputs a comprehensive annotation for splice site variants using hg19 as a reference.
This output includes the following:
  - GERP++2 score
  - Exac03 score
  - 1000 Genomes score
  - MaxEntScan wild-type and variant score
  - MaxEntScan % difference between WT and Mut
  - list of lost conserved domains (rpsblast of Cdd)
  - list of lost pdb sequences (blastp of pdb)
  
  

<a name="installation"/>

## II. Installation Instructions

To install SVA locally follow these steps:
  - clone project
  - from the main directory of the project run maven using `mvn clean -U install`
  - run `java -jar target/SVA-1.0-SNAPSHOT-jar-with-dependencies.jar -help` to see command line options

Required Databases:
  SVA uses different databases and software to provide its results. Below are instructions on how to download and install
  the required databases and software.
  

### Annovar

The main package download for Annovar can be found [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/). Once Annovar software has been downloaded, it can be used to download 3 other required databases. To
download the additional databases use the following command three times:

`./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <name of db> <output_folder(usually humandb/)>`

The abbreviations for the databases used by SVA are as follows:
  - gerp++gt2
  - exac03
  - 1000g2015aug
  
### hg19 Reference Genome

SVA requires the hg19 genome separated by chromosome. Instructions on how to download this are given on UCSC website [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/). The files must be uncompressed as explained in the previous
link.

### UCSC RefSeq file

The required RefSeq file comes from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables). Included in the main directory of the github directory, the user does not need to download their own RefSeq file.
  

<a name="instruct"/>

## III. Usage Instructions and Examples

  
Required command line arguments:

<a name="license"/>

## IV. License

<a name="funding"/>

## V. Funding and Acknowledgements

<a name="contact"/>

## VI. Contact


