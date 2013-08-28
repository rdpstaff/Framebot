This folder contains the protein or nucleotide reference sequences for a few genes. The name convention is <genename>_<prot|nucl>_ref[_primername].fasta for fasta file and <genename>_<prot|nucl>_ref[_primername]_g<number>f<number>e<number>.index for index file. 
-g: gap-open-penalty
-f: frameshift-penalty
-e: gap-ext-penalty

Example:
nifh_prot_ref.fasta: the near full-length nifH protein sequences
nifh_prot_ref_polyprimers.fasta: the protein sequences flanked by the Poly nifH primers (Poly et al. 2001. Appl. Environ. Microbiol. 67:2255-2262.)
nifh_nucl_ref_polyprimers.fasta: the nucleotide sequences flanked by the Poly nifH primers.o

nifh_nucl_ref_polyprimers_g13f10e4.index : the index file built using the nifh_nucl_ref_polyprimers.fasta and paramaters -g -13 -f -10 and -e -4. 


## lineage for the taxonomic classification 
For nifH gene, nifh_prot_ref_lineage.txt contains the description and lineage (based on GenBank annotation) for sequences in file nifh_prot_ref.fasta. This nifh_prot_ref_lineage.txt can be used by the command GetFrameBotStatMain to gather taxonomic classification. See README file in this package.
