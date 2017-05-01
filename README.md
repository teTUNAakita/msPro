<!--
ctrl + shift + M: show preview
-->
# msPro

## Introduction

msPro is a simulation program to generate patterns of single nucleotide polymorphism (SNP) data in a region that experienced homologous recombination with different species. The software was developed by modifying the commonly used Hudson’s ms simulator (Hudson, 1990). msPro is flexible so that it can consider the coalescent process of a particular species (population), around which there are a number of different species and the focal species potentially undergo recombination with these species. msPro also incorporates any change of the population size simultaneously. The software runs coalescent simulations conditional on the joint distribution of divergence and successfully integrated tract length.

## Download and compilation
You can download the most up to date version of msPro via GitHub:
```
git clone https://github.com/TetsuyaAkita/msPro
```
The source code of the program is written in C and this program is intended to be run on UNIX, or a UNIX-like operating systems, such as Linux or Mac OS X. After downloading, change the directory by typing: cd msPro; and then, compile the program by typing:
```
gcc -O3 -Wall -Wextra -o msC msC.c streecC.c rand1.c -lm -w
```

## Running msPro


## Output


- Number of species is limited to ten.
