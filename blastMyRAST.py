from Bio import SeqIO
from Bio import Bio.Blast.Applications

recs = [r for r in SeqIO.parse(myRast, 'genbank')

for r in recs:
    for feature in r.features:
        if feature.type == "CDS":
            print(feature.qualifiers['translation'])

 
