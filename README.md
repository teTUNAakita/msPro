<!--
ctrl + shift + M: show preview
-->
# msPro

## Introduction
**msPro** is a simulation program to generate patterns of single nucleotide polymorphism (SNP) data in a region that have experienced homologous recombination with different species. The software was developed by modifying the commonly used Hudson’s **ms** simulator (Hudson, 1990). **msPro** is flexible so that it can consider the coalescent process of a particular species (population), around which there are a number of different species and the focal species potentially undergoes recombination with these species. **msPro** also incorporates population size changes, population structures, and admixture events simultaneously as **ms** does. The software runs coalescent simulations conditional on the joint probability distribution of divergence and successfully integrated tract length of foreigh DNA. A scientific paper describing **msPro** in detail will be published elsewhere, entitled "The coalescent for prokaryotes with homologous recombination from external source" by Tetsuya Akita, Shohei Takuno, and Hideki Innan.

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
*nsam* is the sample size. *nreps* is the number of independent samples to generate. 2NμL after the ‘-t’ switch is the population mutation parameter per region, where N is the current population size, L is the length of focal region, and μ is the mutation rate per site. ***after the ‘-r’ switch is the crossing-over rate per site and this should be removed?*** Intra-species gene conversion is assumed to initiate at any position at rate *g* per site per generation. Tract length of gene conversion is assumed to follow a geometric distribution with mean λ (Wiuf and Hein, 2000). *2Ng(L-1)* after the `-c` switch is the population gene conversion rate (within species), and λ is the mean conversion tract length. By contrast to intra-species gene conversion, the treatment of inter-species gene conversion is a little tricky. Consider backward argument: the backward initiation is occurred at rate *h* per site per generation (i.e. *h* is the rate at which the lineage experiences a recombination event from external source that is initiated at a site). *2Nh(L-1)* after the `-b` switch is the population gene conversion rate (between species). "dist.txt" specifies the joint probability distribution of divergence and tract length that is successfully integrated, as explained in the next.   

## Format of external resource



## Output

## References
- Hudson (2002) (1990?)

- Wiuf and Hein (2000)

## Memo
- Number of species is limited to ten.

- It's good to add some files of the joint probability distribution of divergence and tract length as examples.

- Test to see if **msPro** runs in Linux.

- Hudson (1990) in Introduction? Or we should cite Hudson (2002) Bioinformatics.  It's better to add the URL to **ms**'s manual.
