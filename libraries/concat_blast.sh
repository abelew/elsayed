#!/bin/env bash
## This script assumes input from make_annotated_genes.pl
## The esmer.fasta in this instance is thelist of 'gene'
## annotations from the Tritryp db gff files...
cp esmer.fasta esmer-non.fasta
cp esmer.fasta esmer-unas.fasta
cp esmer.fasta all.fasta

cp non.fasta non-unas.fasta

cat non.fasta >> esmer-non.fasta
cat unas.fasta >> esmer-unas.fasta
cat unas.fasta >> non-unas.fasta
cat non.fasta >> all.fasta
cat unas.fasta >> all.fasta

## Resulting files:
#esmer.fasta  ## just esmeraldo reads
#non.fasta    ## nonesmeraldo reads
#unas.fasta   ## unassigned contigs
#esmer-non.fasta  ## Both haplotypes
#esmer-unas.fasta ## esmer+unassigned
#non-unas.fasta   ## nonesmer+unassigned
#all.fasta        ## all 3 combined


formatdb -i esmer.fasta -p F -o T -n esmer -s
formatdb -i non.fasta -p F -o T -n non -s
formatdb -i unas.fasta -p F -o T -n unas -s
##  The 3 possible pairs, esmer+non, esmer+unas, non+unas
formatdb -i esmer-non.fasta -p F -o T -n esmer-non -s
formatdb -i esmer-unas.fasta -p F -o T -n esmer-unas -s
formatdb -i non-unas.fasta -p F -o T -n non-unas -s
## And the union of everything
formatdb -i all.fasta -p F -o T -n all -s
