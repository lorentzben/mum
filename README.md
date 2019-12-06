mum
-------------------------------------------------
Performing 16s microbial analysis is made easier with the invention of tools like Mothur, however inputing commands one-by-one is tiring and batch scripts are not generalizable. This script intends to take a directory of raw fastq files and process them following Kelly Lab's general protocol, or with custom input.


## Prerequisities
* Linux
* Mothur v.1.42.2 (installed and in path)
* Python 3 (installed and in path)
* pip3 v.19.2.1
* mothur_py v.0.4.0(pip3 install mothur_py)
* design.txt (User Generated file)
* oligos.txt (User generated file)

## Install

```shell
$ git clone git@github.com:lorentzben/mum.git
```
TODO fill out what else will be in here
After cloning a folder called mum will be created. Inside will be this README and _____
## Help
```shell
$ python3 mum.py -h 

usage: mum.py [-h] -n JOB_NAME [-c] [-l MAX_LEN] [-p PRE_CLUST]
              [-s SUB_SAMPLE] [-v]

Perform Analysis of 16s Microbial Data using Mothur

optional arguments:
  -h, --help     show this help message and exit
  -n JOB_NAME    name for analysis, will be filename for resultant files
  -c             flag to use custom parameters, ignore for kelly lab
                 parameters
  -l MAX_LEN     determines the longest value permitted for sequences
  -p PRE_CLUST   pre cluster value, higher is more stringent
  -s SUB_SAMPLE  sub sample value, only applicable for custom runs
  -v, --version  show program's version number and exit
```

mum can be calles with only a name given, or with custom running parameters
## Running 
```shell
$ chmod +x mum.py
$ python3 mum.py -n NAME
OR
$ python3 mum.py -n NAME -c -l "###" -p "###" -s "###"
```
Ensure that the files design.txt, oligos.txt as well as the directory of either fastq files or nested files from illumina are all in the current directory. The results of analysis will be placed ____ 

TODO fill in where this location is

Running
```shell
$ python3 cleanup.py
```
This command will remove some files created by mum, however the user will still need to run:
```shell
$ rm NAME*
```
## Output
TODO fill in names

Inside of the directory ______there will be _______

## Current Files
* mum.py
* mothur_batch.py
* MothurReader.py (Diagnositc File) 
* setup.py (Diagnostic File)
* cleanup.py (Diagnostic File)
* README.md (This File)
* design.txt (example file)
* oligos.txt (example file)


## Version
* Version 1.0

## Author
* Ben Lorentz

## Future Plans
* Implement verbose logging at the info and debugging scale
* Implememt unit testing using the unittest package
* Move final files into a output directory 
