# Mum
This is a tool to aid in the analysis of 16s data using the tool Mothur

This tool comes with a few options:
1. Custom Filename for the Mothur Analysis (Required)
2. Option to follow Kelly Lab MiSeq protocol for Mothur or Option to input custom parameters to be more or less stringent. 

Required:
* Linux
* Mothur (installed and in path)
* Python 3 (installed and in path)
* mothur_py (pip3 install mothur_py)
* design.txt (User Generated file)
* oligos.txt (User generated file)

The current files are:
* mum.py
* mothur_batch.py
* MothurReader.py (Diagnositc File) 
* setup.py (Diagnostic File)
* cleanup.py (Diagnostic File)
* README.md (This File)
* design.txt (example file)
* oligos.txt (example file)

Future work:
* Implement verbose logging at the info and debugging scale
* Implememt unit testing using the unittest package
* Move final files into a output directory 
