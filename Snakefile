import os
import glob

configfile: "config.yaml"

PHAGE = [os.path.basename(f)[0] for f in glob.glob('*.gbk')]

#FINAL = expand("mafft/{sample}_aligned.fasta", sample=SAMPLE)

rule all:
    input: FINAL

rule convert_genbank_to_fasta:
    input: '{myrast}.gbk'
    output: '{myrast}.fasta'
    script: 'blastMyRAST.py'

rule psiblast:
    input: '{myrast}.fasta'
    output: 'psiblast/{phage}_psiblast.tab'
    shell: 'psiblast -i {input} -db {config[blastdb]} -o {output}'

rule transterm_terminator:
    input: '{myrast}.fasta'
    output: 'transterm/{transterm}.tt'
    shell: 'transterm -p expterm.dat {input} {config[annotation]}.coords > {output}


rule compare_annotations:
    input: ''
    output: 'mafft/{sample}_aligned.fasta'
    shell: 'mafft {input} > {output}'

rule aragorn:
    input: '{myrast}.fasta'
    output: 'aragorn_tRNA/{myrast}_aragorn.fasta'
    shell: 'aragorn -t -o {output} {input}

# Potentially add tRNAscan-SE if standalone version available

# TODO: run trial for checking output file extension
rule phobius:
    input: '{myrast}.fasta'
    output: 'transmembrane/phobius/{myrast}.txt'
    shell: 'phobius.pl {input}'

# TODO: receive the program from CBS Denmark
rule tmhmm:
    input: '{myrast}.fasta'
    output: 'transmembrane

