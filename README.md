# **PolyGEN**

Design automation of polycistronic tRNA-based genes containing custom RNAs for assembly in BsaI-driven Golden Gate experiments. The entire backbone of PolyGEN is based heavily on
[iBioCAD](https://ibiocad.igb.illinois.edu/) by HamediRad et al. (2019), for which the code is openly available [here](https://github.com/scottweis1/iBioCAD).

The code takes as input an array of custom RNAs and will compute the finished PTG together with the necessary oligomers to produce all parts from a plasmid containing a gRNA-tRNA template. Currently,
the produced PTGs can include sgRNAs, pegRNAs and other small RNAs. By default, PolyGEN will use the following parameters, which can be varied

- primer melting temperature between 52 and 72 °C if possible
- primer annealing length on template between 12 and 30 nt
- 'tgcc' and 'gttt' as BsaI restriction overlaps with the plasmid
- no additional BsaI restriction sites in the plasmid

To calculate the primer melting temperatures, PolyGEN uses the same method and parameters as Benchling: [SantaLucia (1998)](https://www.pnas.org/content/95/4/1460).

PolyGEN can be accessed as a webapp:

CRISPR and small RNA technologies 
