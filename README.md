SpliceSiteVariantAnalyzer (SSVA)

#### Table of Contents  
[I. Introduction](#introduction)  
[II. Installation Instructions](#installation)  
[III. Usage Instructions and Examples](#instruct)        
[IV. Contact](#contact)   


<a name="introduction"/>   

## I. Introduction

SpliceSiteVariantAnalyzer (SSVA) is a splice site variant diagnosis tool which outputs a comprehensive annotation for splice site variants from a vcf file. The software is compatible with both hg19 or hg38 as the reference.

This output includes the following:
  - GERP++2 score
  - Exac03 score
  - 1000 Genomes score
  - dbSCSNV score
  - MaxEntScan wild-type and variant score
  - MaxEntScan % difference between WT and Mut
  - list of lost conserved domains (rpsblast of Cdd)  
  

<a name="installation"/>

## II. Installation Instructions

To install SSVA locally follow these steps:
  - clone project
  - from the main directory of the project run maven using `mvn clean -U install`
  - run `java -jar target/SSVA-1.0-jar-with-dependencies.jar -help` to see command line options

Required Databases:
  SSVA uses different databases and software to provide its results. Below are instructions on how to download and install
  the required databases and software.
  

### Annovar

The main package download for Annovar can be found [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/). Once Annovar software has been downloaded, it can be used to download other required databases. To
download the additional databases use the following command for each database:

`./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <name of db> <output_folder(usually humandb/)>`

The abbreviations for the databases (\<name of db\>) used by SSVA are as follows:
  - refGene
  - dbscsnv11
  - gerp++gt2 (hg19 only)
  - exac03
  - 1000g2015aug
  
### Reference Genome (hg19/hg38)

SSVA requires either the hg19 or hg38 genome separated by chromosome. Instructions on how to download this are given on UCSC website [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/). The files must be uncompressed as explained in the previous
link.

### UCSC RefSeq file

Both the hg38 and hg19 RefSeq files are included in the main directory of the github repository, therefore the user does not need to download their own RefSeq file. These files comes from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables).

### MaxEntScan scripts

The MaxEntScan scripts are already downloaded in the git repository under the `MES` and `splicemodels` directories, therefore the user does not need to download the MES scripts. These scripts were downloaded from https://github.com/razZ0r/maxentscan.git.
  

<a name="instruct"/>

## III. Usage Instructions and Examples

Required command line arguments:
- -i input file:  The vcf file with all variants
- -o output directory:  The directory where the output files will be saved (the directory must exist already)
- -a annovar: The path to the annovar directory with perl scripts 
- -d humandb:  The path to the directory containing the databases that uses (must have refGene, dbscsnv11, gerp++gt2 (hg19 only), exac03, and 1000g2015aug as explained above)
- -g genome:  The path to the directory that contains the UCSC reference genome by chromosome downloaded (hg19/hg38)

Optional command line arguments:
- -b build version:  Valid options are 'hg19'(default) or 'hg38'
- -r RefSeqFile:  This is the path to the file that contains the UCSC table viewer RefSeq data. (default: in this git repository)
- -m MaxEntScan path:  The path to the MaxEntScan perl scripts directory (default: ./MES/)
- -e evalue:  This is the evalue cut off for the rpsblast of the Cdd database. (default: .005)
- -s samtools:  This is the path to the Samtools executable. (default: in path variable)
- -p rpsblast: This toggles whether to do rpsblast step or not. (default: false)

An example for running the program is as follows:

```
java -jar ./target/SSVA-1.0-jar-with-dependencies.jar \
        -i chromosome22.vcf \
        -o SVA_output/ \
        -a ~/software/annovar \
        -d ~/humandb/ \
        -g ~/hg19Ref/ \
```


<a name="contact"/>

## IV. Contact

For questions, comments, concerns, feature requests, suggestions, etc., please
contact:

Pery Ridge, Ph.D. -- perry.ridge@byu.edu

Note: For usage questions, please consult section `III. Usage Instructions and
Examples' first.
