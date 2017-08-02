import os
import glob

configfile: "config.yaml"

PHAGE = [os.path.basename(f)[0] for f in glob.glob('*.gbk')]

#FINAL = expand("mafft/{sample}_aligned.fasta", sample=SAMPLE)

rule all:
    input: FINAL

rule genbank_to_aa_fasta:
    input: '{myrast}.gbk'
    output: '{myrast}.fasta'
    script: 'blastMyRAST.py'

#Need to write script for this, or add function to script above
rule genbank_to_genome_fasta:
    input: '{myrast}.gbk'
    output: '{myrast}_complete_genome.gbk'

# Homolog search with psiblast
rule psiblast:
    input: '{myrast}.fasta'
    output: 'psiblast/{phage}_psiblast.tab'
    shell: 'psiblast -i {input} -db {config[blastdb]} -o {output}'

# Promoter prediction: need to download neural network software

# Search for rho-dependent terminators
rule transterm_terminator:
    input: '{myrast}.fasta'
    output: 'transterm/{transterm}.tt'
    shell: 'transterm -p expterm.dat {input} {config[annotation]}.coords > {output}

#rule compare_annotations:
#    input: ''
#    output: 'mafft/{sample}_aligned.fasta'
#    shell: 'mafft {input} > {output}'

# tRNA screening

rule aragorn:
    input: '{myrast_dna}.fasta'
    output: 'tRNA_screening/{myrast}_aragorn.fasta'
    shell: 'aragorn -t -o {output} {input}

rule tRNAscan:
    input: '{myrast_dna}.fasta
    output: 'tRNA_screening/{myrast_dna}.fasta'
    shell: tRNAscan-SE -qQ -Y -o# -m# -f# -l# -c tRNAscan-SE.conf -s# (fasta file)

# TODO: run trial with toy data for checking output file extension
rule phobius:
    input: '{myrast}.fasta'
    output: 'transmembrane/phobius/{myrast}.txt'
    shell: 'phobius.pl {input}'

# TODO: receive the program from CBS Denmark and check usage
rule tmhmm:
    input: '{myrast}.fasta'
    output: 'transmembrane/TMHMM/{myrast}.txt'
    shell: 'tmhmm {input} > {output}'

# TODO: CLUSTALW
rule clustalw:
    input:
    output:
    shell:
    

