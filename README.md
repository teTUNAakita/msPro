<!--
ctrl + shift + M: show preview
-->
# msPro

## Introduction
**msPro** is a simulation program to generate patterns of single nucleotide polymorphism (SNP) data in a region that have experienced homologous recombination with different species. The software was developed by modifying the commonly used Hudson’s **ms** simulator (Hudson, 2002; URL: http://home.uchicago.edu/rhudson1/source/mksamples.html). **msPro** is flexible so that it can consider the coalescent process of a particular species (population), around which there are a number of different species and the focal species potentially undergoes recombination with these species. **msPro** also incorporates population size changes, population structures, and admixture events simultaneously as **ms** does. The software runs coalescent simulations conditional on the joint probability distribution of divergence and successfully integrated tract length of foreigh DNA. A scientific paper describing **msPro** in detail will be published elsewhere, entitled "The coalescent for prokaryotes with homologous recombination from external source" by Tetsuya Akita, Shohei Takuno, and Hideki Innan.

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
./msPro nsam nreps L -t 2NμL -c 2Ng(L-1) λ -b 2Nh(L-1) dist.txt
```
*nsam* is the sample size. *nreps* is the number of replicates to generate. *L* is the length of a focal region. 2*NμL* after the ‘-t’ switch is the population mutation parameter per region, where *N* is the current population size, and μ is the mutation rate per site per generation. Intra-species gene conversion is assumed to initiate at any position at rate *g* per site per generation. Tract length of gene conversion is assumed to follow a geometric distribution with mean *λ* (Wiuf and Hein, 2000). 2*Ng(L-1)* after the `-c` switch is the population gene conversion rate (within species), and λ is the mean conversion tract length. Treatment of inter-species gene conversion is based on backward argument: the backward initiation is occurred at rate *h* per site per generation (i.e. *h* is the rate at which the lineage experiences a recombination event from external source that is initiated at a site). 2*Nh(L-1)* after the `-b` switch is the population gene conversion rate (between species). "dist.txt" specifies the name of a file containing joint probability distribution of divergence and tract length that is successfully integrated (see Fig. XX in our paper). This prior distribution is necessary for running **msPro** and is located in "msPro" directory.

Here is an example of a command line:
```
msPro 5 2 50 -t 1.0 -c 1 10 -b 0.1 dist.txt
```
In this case, the program generate two data sets of five DNA sequences. The parameters were set as 2*NμL* = 1.0, 2*Ng(L-1)* = 1.0, *λ* = 10, 2*Nh(L-1)* = 0.1, and the prior distribution of divergence and tract length is specified in "dist.txt", as explained in the following section.

## Format of a prior distribution of divergene and tract length of foreign DNAs
An example of the format of a input file (specified as "dist.txt" in this case) is as follows:
```
10	0.3	0.05
50	0.2	0.5
100	0.1	0.3
1000	0.05	0.15
```

This file indicates that there are four cases of integrated foreign DNAs. The first column is the tract length of integrated DNAs (*x*). The second column is the average nucleotide divergence between the simulated species and one of the species as sourses of foreign DNAs (*d*)；*d* works as a species name. The third column is the frequency given *x* and *d* (*f*). *x* > 0 and 0 < *d* < 1 must be satisfied. The frequency is conditioned such that its summation becomes unity. The number of species around the focal population (i.e., *d*) is limited to ten.  

## Output

The output format is as follows :

```
./msPro 5 2 50 -t 1 -c 1 10 -b 0.1 dist.txt
16723 65395 1320

//
segsites: 3, nch: 3
00000000010000000000010000000000010000000000000000
00000000000001000000000000000000000000000000000000
00000000000001000000000000000000000000000000000000
00000000000001000000000000000000000000000000000000
00000000000001000100000000000000000000000000000000


//
segsites: 2, nch: 2
01000000000000001010011001001000000000011000000000
00001000100000000000000000000000000100000000000000
01000000000000001010011001001000000000011000000000
00001000100000000000000000000000000100000000000000
01000000000000001010011001001000000000011000000000
```
The first line of the output is the command line.
The second line shows the random number seeds.
The output contains two replicates in this case.
Each replicate starts with "segsites: X, nch: Y", where X is the number of mutations events and Y is the number of inter-species recombination evens.
Following this lines, you find five chromosomes in each replicate, and 0 and 1 represent different nucleotides.
It should be noted that **msPro** allows back mutation events, and therefore the allele "1" is not necessarily a derived allele???


## References
- Hudson, R R., (2002) Generating samples under a Wright-Fisher neutral model. Bioinformatics 18:337-338,

- Wiuf, C., and Hein, J., (2000) The coalescent with gene conversion. Genetics 155:451-462

## Memo
<!--
- Number of species is limited to ten.

- It's good to add some files of the joint probability distribution of divergence and tract length as examples.

- Test to see if **msPro** runs in Linux.

- Hudson (1990) in Introduction? Or we should cite Hudson (2002) Bioinformatics.  It's better to add the URL to **ms**'s manual.
-->
- Test to see if **msPro** runs in Linux.

- Think about we should submit our paper to bioRxiv **TA: I agree**

- Add explanation of important options.  
