<!--
ctrl + shift + M: show preview
-->
# msPro

## Introduction
**msPro** is a simulation program to generate patterns of single nucleotide polymorphism (SNP) data in a region that experienced homologous recombination with different species. The software was developed by modifying the commonly used Hudson’s **ms** simulator (Hudson, 1990). **msPro** is flexible so that it can consider the coalescent process of a particular species (population), around which there are a number of different species and the focal species potentially undergo recombination with these species. **msPro** also incorporates any change of the population size simultaneously. The software runs coalescent simulations conditional on the joint distribution of divergence and successfully integrated tract length. A scientific paper describing **msPro** in detail will be published elsewhere, entitled "The coalescent for prokaryotes with homologous recombination from external source" by Tetsuya Akita, Shohei Takuno, and Hideki Innan.

## Download and compilation
You can download the most up to date version of msPro via GitHub:
```
git clone https://github.com/TetsuyaAkita/msPro
```
The source code of the program is written in C and this program is intended to be run on UNIX, or a UNIX-like operating systems, such as Linux or Mac OS X. After downloading, change the directory by: `cd msPro`; and then, compile the program:
```
gcc -O3 -Wall -Wextra -o msPro msPro.c streecPro.c rand1.c -lm -w
```

## Running msPro
The following command line shows the simplest usage of msPro:
```
./msC nsam nreps -t 2NμL -r 0.0 L -c 2Ng(L-1) λ -b 2Nh(L-1) dist.txt
```
*nsam* is the sample size. *nreps* is the number of independent samples to generate. *2NμL* after the `-t` switch is the population mutation parameter per region, where N is the current population size, *L* is the length of focal region, and *μ* is the mutation rate per site. ***0.0 after the ‘-r’ switch is the crossing-over rate per site and this should be removed.*** Gene conversion is assumed to initiate between a specified pair of sites in a given chromosome with probability *g*, and the conversion tract length is assumed to follow geometric distribution with mean λ (Wiuf and Hein, 2000). *2Ng(L-1)* after the `-c` switch is the population gene conversion rate.
λ is the mean conversion tract length.   
## Output

## Memo
- Number of species is limited to ten.
