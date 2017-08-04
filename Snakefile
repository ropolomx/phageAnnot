import os

# Globals 
configfile: 'config.yaml'

# Path to the directory with Genbank files
SAMPLEDIR = config['genbank']

# Names of bacteriophages
PHAGES = [os.path.splitext(os.path.basename(f))[0] for f in os.listdir(SAMPLEDIR) if f.endswith('.gbk')] 

# Names of BLAST and ARAGORN targets

BLAST = expand("psiblast/{phage}_psiblast.tab", phage=PHAGES)
ARAGORN = expand("tRNA_screening/{phage}_aragorn.fasta", phage=PHAGES)
TMHMM = expand("transmembrane/TMHMM/{phage}_tmhmm.txt", phage=PHAGES)

rule all:
    input: BLAST, ARAGORN, TMHMM

# TODO: Make quality checks, pre-processing, assembly, and annotation optional

# Generate FASTQC reports

#rule fastqc_analysis:
#    input: 'config["samples"]/{reads}.fastq'
#    output: 'quality/fastqc_{reads}/{reads}_fastqc.html', dir='quality/fastqc_{reads}'
#    shell: 'fastqc {input} -o {output.dir}'

# Aggregate FASTQC reports with MultiQC
# TODO: Fix the input names. Maybe with expand...

#rule multiqc:
#    input: 'quality/fastqc_*/*_fastqc.html'
#    output: 'quality/multiqc/multiqc_report.html', dir='quality/multiqc'
#    shell: 'multiqc -o {output.dir}'

# Assemble with SPAdes

#rule spades:
#    input: '{sample}.fastq'
#    output: 'spades/{sample}_spades/contigs.fasta', dir='spades/{sample}_spades'
#    shell: 'spades.py -s {input} -o {output.dir}'

# Annotate assemblies with phage

#rule phage:
#    input: 'spades/{sample}_spades/contigs.fasta'
#    output: 'phage/{sample}.gbk'
#    shell: 'phage'

#TODO: fix output of conversion script

rule genbank_to_fastas:
    input: 'Genbank/{phage}.gbk'
    output: aa='fastas/{phage}_amino.fasta', nt='fastas/{phage}_genome.fasta'
    shell: './annotationToFASTA.py -a {output.aa} -g {output.nt} {input}'

# Homolog search with psiblast

rule psiblast:
    input: 'fastas/{phage}_amino.fasta'
    output: 'psiblast/{phage}_psiblast.tab'
    shell: 'psiblast -query {input} -db {config[blastdb]} -outfmt 6 -out {output}'

# Promoter prediction: need to download neural network software

# Search for rho-dependent terminators

#rule transterm_terminator:
#    input: '{phage}_aa.fasta'
#    output: 'transterm/{transterm}.tt'
#    shell: 'transterm -p expterm.dat {input} {config[annotation]}.coords > {output}'

# tRNA screening

rule aragorn:
    input: 'fastas/{phage}_genome.fasta'
    output: 'tRNA_screening/{phage}_aragorn.fasta'
    shell: 'aragorn -fo -o {output} {input}'

#rule tRNAscan:
#    input: '{phage_dna}.fasta'
#    output: 'tRNA_screening/{phage_dna}.fasta'
#    shell: 'tRNAscan-SE -qQ -Y -o# -m# -f# -l# -c tRNAscan-SE.conf -s# (fasta file)'

# Transmembrane domains

#TODO: run trial with toy data for checking output file extension
#rule phobius:
#    input: '{phage}.fasta'
#    output: 'transmembrane/phobius/{phage}.txt'
#    shell: 'phobius.pl {input}'

rule tmhmm:
    input: 'fastas/{phage}_amino.fasta'
    output: 'transmembrane/TMHMM/{phage}_tmhmm.txt'
    shell: 'tmhmm {input} > {output}'

#rule final_multiqc_report:
#    input: ''
#    output: ''
#    shell: 'multiqc .'
