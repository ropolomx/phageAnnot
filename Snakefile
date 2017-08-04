import os

configfile: 'config.yaml'

SAMPLEDIR = config['genbank']

PHAGES = [os.path.splitext(os.path.basename(f))[0] for f in os.listdir(SAMPLEDIR) if f.endswith('.gbk')] 

BLAST = expand("psiblast/{phage}_psiblast.tab", phage=PHAGES)

rule all:
#    input: "report.html"
    input: BLAST

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

# Annotate assemblies with myRAST

#rule myrast:
#    input: 'spades/{sample}_spades/contigs.fasta'
#    output: 'myRAST/{sample}.gbk'
#    shell: 'myrast'

#TODO: fix output of conversion script

rule genbank_to_fastas:
    input: 'Genbank/{myrast}.gbk'
    output: aa='fastas/{myrast}_aa.fasta', nt='fastas/{myrast}_genome.fasta'
    shell: 'annotationToFASTA.py {input}'

# Homolog search with psiblast
rule psiblast:
    input: '{phage}_aa.fasta'
    output: 'psiblast/{phage}_psiblast.tab'
    shell: 'psiblast -i {input} -db {config[blastdb]} -outfmt 6 -o {output}'

# Promoter prediction: need to download neural network software

# Search for rho-dependent terminators

#rule transterm_terminator:
#    input: '{myrast}_aa.fasta'
#    output: 'transterm/{transterm}.tt'
#    shell: 'transterm -p expterm.dat {input} {config[annotation]}.coords > {output}'

# tRNA screening

#rule aragorn:
#    input: '{myrast_dna}.fasta'
#    output: 'tRNA_screening/{myrast}_aragorn.fasta'
#    shell: 'aragorn -t -o {output} {input}'

#rule tRNAscan:
#    input: '{myrast_dna}.fasta'
#    output: 'tRNA_screening/{myrast_dna}.fasta'
#    shell: 'tRNAscan-SE -qQ -Y -o# -m# -f# -l# -c tRNAscan-SE.conf -s# (fasta file)'

# Transmembrane domains

# TODO: run trial with toy data for checking output file extension
#rule phobius:
#    input: '{myrast}.fasta'
#    output: 'transmembrane/phobius/{myrast}.txt'
#    shell: 'phobius.pl {input}'

#rule tmhmm:
#    input: '{myrast}.fasta'
#    output: 'transmembrane/TMHMM/{myrast}.txt'
#    shell: 'tmhmm {input} > {output}'

#rule final_multiqc_report:
#    input: ''
#    output: ''
#    shell: 'multiqc .'
