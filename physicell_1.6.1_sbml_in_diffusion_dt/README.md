# PhysiCell_1.6.1_sbml


## Overview: 
Demonstrate using PhysiCell togetheer with an intracellular (metabolic) model (represented as a SBML file) which is solved using libRoadRunner.

## Dependencies
* modern g++ compiler (rf. Quickstart)
* Python 3 (Anaconda distribution recommended)

## Instructions
```
$ cd beta
$ python setup_libroadrunner.py
```
Follow instructions from that Python script, including editing the PhysiCell Makefile. Then:
```
$ make
$ ftest
```
