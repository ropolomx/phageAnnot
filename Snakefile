import os
import glob

MYRAST = [os.path.basename(f)[0] for f in glob.glob('*.fastq')]

FINAL = expand("mafft/{sample}_aligned.fasta", sample=SAMPLE)

rule all:
    input: FINAL

rule psiblast:
    input: '{myrast}.gbk'
    output: 'spades/{sample}_spades/contigs.fasta', dir='spades/{sample}_spades'
    shell: 'psiblast -i {input} -db -o {output.dir}'

rule compare_annotations:
     input: ''
     output: 'mafft/{sample}_aligned.fasta'
     shell: 'mafft {input} > {output}'
