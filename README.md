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

SpliceVariantAnalyzer (SVA) is a splice site variant diagnosis tool which outputs a comprehensive annotation for splice site variants.
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

The main package download for Annovar can be found [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/). Once Annovar software has been downloaded, it can be used to download 4 other required databases.
    
    
<a name="instruct"/>

## III. Usage Instructions and Examples

  
Required command line arguments:

<a name="license"/>

## IV. License

<a name="funding"/>

## V. Funding and Acknowledgements

<a name="contact"/>

## VI. Contact


