# Riboflavin-gene annotation and pangenome tree

## Dependencies

- [Install eggnog mapper](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Installation)
- [Install scarap](https://github.com/swittouck/scarap#installation)
- [Install iqtree](http://www.iqtree.org/doc/Quickstart) 

## Functional annotation of the isolates with eggnog

To map the genomes against the eggnog database, the following line can be run:
```./src/riboflavin_genes/run_eggnog.sh <list-of-assemblies>```
Where <list-of-assemblies> is a text file containing the path to a genome of interest per line. The command will generate three files per genome:

- <genome>.emapper.annotations
- <genome>.emapper.hits
- <genome>.emapper.seed_orthologs
More information on these files [can be found here](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#output-files)

## Build a phylogenetic tree of the assemblies

Running the command 
```./src/riboflavin_genes/build_tree.sh <faa-file-of-genome>```
On a file containing one .faa file per line will construct a phylogenetic tree of the genomes in the file  

## Analysis of the results

The annotations from the eggnog database are then collected and processed using two quarto script (src/list_riboflav.qmd and src/list_fleet.qmd).
