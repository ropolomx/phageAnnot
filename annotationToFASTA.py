#! /usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse


def arguments():
    parser = argparse.ArgumentParser(description='Convert GenBank CDS to AA and genome to FASTA')
    parser.add_argument('-a','--amino', required=True, help='Name of aminoacid FASTA file output')
    parser.add_argument('-g','--genome', required=True, help='Name of nucleotide genome FASTA file output')
    parser.add_argument('annotation', help='GenBank-formatted annotation file from myRAST or Prokka')
    return parser.parse_args()

# Parse Genbank file

def readGenbank(myRast):

    """
    Read and Parse Seq Records from Genbank file
    """

    myRASTGB = [r for r in SeqIO.parse(myRast, 'genbank')]

    return myRASTGB

def parseGenbankCDS(recs):

    """
    Parse sequence records headers, descriptions and translated CDS sequences from
    the Genbank file
    """

    headers = []

    descriptions = []

    translations = []

    for r in recs:
        for feature in r.features:
            if feature.type == "CDS":
                headers.append('_'.join([r.name, feature.qualifiers['locus_tag'][0]]))
                descriptions.append(feature.qualifiers['product'][0])
                translations.append(feature.qualifiers['translation'][0])

    fastaTuple = zip(headers, descriptions, translations)

    return fastaTuple

def writeAminoFasta(fastaInfo, outputAminoFile):

    """
    List comprehension to generate new sequence records from the information
    generated in the tuple above. Returns writing sequence records into an aminoacid
    FASTA file
    """
    newRecs = [SeqRecord(Seq(f[2], IUPAC.protein),
        id=f[0],
        name=f[0],
        description=f[1]) for f in fastaInfo]

    # Write new records as FASTA
    return SeqIO.write(newRecs, outputAminoFile,'fasta')

def writeGenomeFasta(recs, outputGenomeFile):

    """
    Writes genome sequence to a nucleotide FASTA file
    """
# Keep an eye on Genbank records with multiple genomes

    return SeqIO.write(recs, outputGenomeFile, 'fasta')


def main():

    args = arguments()

    genbankRecs = readGenbank(args.annotation)

    fastaElements = parseGenbankCDS(genbankRecs)

    writeAminoFasta(fastaElements, args.amino)

    writeGenomeFasta(genbankRecs, args.genome)

if __name__ == '__main__':
    main()
