# VirSeqImprover

An integrated pipeline for error-correction, extension, and annotation of viral scaffolds to recover error-free full genomes of viruses and phages by polishing and extending draft viral assemblies.

# Installation

Use conda install to install VirSeqImprover.

```bash
conda install virseqimprover
```

# Run VirSeqImprover

## Input

A viral scaffold/contig (.fasta) and the read sample (.fastq) from which the input scaffold is assembled.

The reads can be paired-end reads or single-end reads. For paired-end reads, two separate read files are needed.

## Usage

## Output

The ouput folders and files will be generated in the output directory folder you specify (the same directory of your input files by default if not specified).

The final output file -- ```pilon_out.fasta``` is the final improved scaffold of the original scaffold.

# Annotation

After applying VirSeqImprover to the viral scaffold, use eggNOG-mapper (http://eggnog-mapper.embl.de/) to annotate viral genes. On the tool website, choose "Metagenomic" as the data kind, Prodigal as the gene prediction method, upload the improved scaffold file, and choose "Viruses - 10239" as the taxonomic scope.

When the job is done, download the ```out.emapper.decorated.gff``` file, rename it as ```decorated.gff```.

Use ```decorated.gff``` and ```pilon_out.fasta``` to generate the annotation sheets by running ```ann.py```.
