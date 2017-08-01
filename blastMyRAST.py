#! /usr/bin/env python3


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Parse Genbank file

myRASTGB = [r for r in SeqIO.parse(myRast, 'genbank')

# Converting CDS regions from Gebank file to FASTA for blastp/psiblast search

headers = []

descriptions = []

translations = []

for r in recs:
    for feature in r.features:
        if feature.type == "CDS":
            headers.append('_'.join([r.name, feature.qualifiers['locus_tag'][0]) 
            descriptions.append(feature.qualifiers['product'][0]])) 
            translations.append(feature.qualifiers['translation'][0])


# Zip headers, descriptions, and translated sequences into a zip object with tuples
# Returns an iterator of tuples

fastaTuple = zip(headers, descriptions, translations) 

# List comprehension to generate new sequence records from the information
# generated in the tuple above

# TODO: fix description

newRecs = [SeqRecord(Seq(f[2], IUPAC.protein), id=f[0], name=f[0], description=f[1]) for f in fastaTuple]

# Write new records as FASTA

SeqIO.write(newRecs, 'newRecs.fasta', 'fasta')
