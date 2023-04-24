import sys
import getopt
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

def gc(dna):
    "This command takes a seq as a string and returns the gc percentage of it."
    for i in dna:
        if i not in "AGCTNagctn":
            return "INVALID SEQUENCE"
    nbases = dna.count('n') + dna.count('N')
    gcpercent = float(dna.count('c') + dna.count('C') + dna.count('g')
                      + dna.count('G')) / (len(dna) - nbases)*100
    return "GC Percentage="+ str(gcpercent) + "%"

#------------------------------------------------------------------------
def transcribe(sequence):
    "This command takes a seq as a string and returns its transcription."
    for i in sequence:
        if i not in "AGCTagct":
            return "INVALID SEQUENCE"
    seq = Seq(sequence)
    return seq.transcribe()

#------------------------------------------------------------------------

def reverse_complement(sequence):
    "This command takes a seq as a string and returns its reverse complement."
    for i in sequence:
        if i not in "AGCTagct":
            return "INVALID SEQUENCE"
    seq = Seq(sequence)
    return seq.reverse_complement()

#------------------------------------------------------------------------

def calc_nbases(sequence):
    "This command takes a seq and calculates its nbases"
    for i in sequence:
        if i not in "AGCTNagctn":
            return "INVALID SEQUENCE"
    return sequence.count('n')+sequence.count('N')

#------------------------------------------------------------------------

def is_valid(sequcence,sequenceType):
    """This command takes a seq and a type (protein, dna, rna) and returns a
    Boolean value of whether it’s a valid type or not"""
    dna = "AGCT"
    rna = "AGCU"
    protein = "ABCDEFGHIKLMNPQRSTVWXYZ"
    sequcence = sequcence.upper()
    sequenceType = sequenceType.lower()
    if sequenceType == "protein":
        for i in sequcence:
            if i not in protein:
                return "Not Valid"
        else:
            return "Valid"
    elif sequenceType == "dna":
        for i in sequcence:
            if i not in dna:
                return "Not Valid"
        else:
            return "Valid"
    elif sequenceType == "rna":
        for i in sequcence:
            if i not in rna:
                return "Not Valid"
        else:
            return "Valid"
    else:
        return "Invalid sequence type"

#------------------------------------------------------------------------

def filter_nbases(sequece):
    "This command takes a seq and returns the Seq after removing n bases"
    newSequence = []
    for i in sequece:
        if i == "N" or i == "n":
            pass
        else:
            newSequence.append(i)
    return "".join(newSequence)

#------------------------------------------------------------------------

def seq_alignment(seq1,seq2, file_name=""):
    """This command takes 2 sequences as input and get all their alignments
       along with the score and the -o is an optional parameter
       if we need the output to be written on a file instead of the screen"""
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    aligns = pairwise2.align.globalxx(seq1, seq2)
    if(len(file_name)==0):
        for i in aligns:
            print(format_alignment(*i))
    else:
       x=file_name.rfind("txt")
       if(x==len(file_name)-3):
          f=open(file_name,"w")
          for k in aligns:
            f.write(format_alignment(*k))
          print("file is created successfully")
       else:
          print("file should be text file")

#------------------------------------------------------------------------

def seq_alignment_files(s1,s2,file_name=''):
    """This command takes 2 fasta files as input, each file contains a single
   sequence. It reads the 2 sequences from files and get all their alignments
   along with the score. The -o is an optional parameter if we need the output
   to be written on a file instead of the screen."""

    record1 = SeqIO.read(s1, "fasta")
    record2 = SeqIO.read(s2, "fasta")
    aligns = pairwise2.align.globalxx(record1.seq, record2.seq)
    if(len(file_name)==0):
        for i in aligns:
            print(format_alignment(*i))
    else:
       x=file_name.rfind("txt")
       if(x==len(file_name)-3):
          f=open(file_name,"w")
          for k in aligns:
            f.write(format_alignment(*k))
          print("file is created successfully")
       else:
          print("file should be text file")

#------------------------------------------------------------------------

def online_alignment(seq , file=""):
    result_handle = NCBIWWW.qblast("blastn", "nt", seq)
    blast_record = NCBIXML.read(result_handle)
    if len(file)==0:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)

    else:
            f = open(file, "w")
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:

                    f.write("****Alignment****\n")
                    f.write('sequence:')
                    f.write(alignment.title)
                    f.write('\n')
                    f.write('length:')
                    f.write(str(alignment.length))
                    f.write('\n')
                    f.write('e value:')
                    f.write(str(hsp.expect))
                    f.write('\n')
                    f.write(hsp.query)
                    f.write('\n')
                    f.write(hsp.match)
                    f.write('\n')
                    f.write(hsp.sbjct)
                    f.write('\n')
                    f.write("---------------------------")
                    f.write('\n')
            f.close()

#------------------------------------------------------------------------

def merge_fasta(*paths,outputPath=''):
    """This command merges fasta files the -o is an optional parameter
       if we need the output to be written on a file instead of the screen"""
    for filename in paths:
        with open(filename) as input_handle:
            sequences = SeqIO.parse(input_handle, "fasta")
            if len(outputPath) != 0:
                with open(outputPath, "a") as output_handle:
                     SeqIO.write(sequences, output_handle, "fasta")
            else:
                for seq_record in SeqIO.parse(input_handle, "fasta"):
                    print(seq_record.description)
                    print(seq_record.seq)

#------------------------------------------------------------------------

def convert_to_fasta(fileName):
    """This command converts the input genbank file with multiple records onto a
    fasta formatted file. The output is to be written in a different output fasta file."""
    if fileName[-3:] == "gbk":
        with open(fileName) as inputFile, open(fileName[:-3]+"fasta", "w") as outputFile:
            sequences = SeqIO.parse(inputFile,"genbank")
            SeqIO.write(sequences,outputFile,"fasta")
            return "File converted successfully"
    else:
        return "Must be a genbank file"

#------------------------------------------------------------------------

def main():
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],  'o:')
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    if args != []:
        if args[0] == 'gc':
            if len(args) == 2:
                print(gc(args[1]))
            else:
                print('\'%s\' command expects one argument'%args[0])
        elif args[0] == "transcribe":
            if len(args) == 2:
                print(transcribe(args[1]))
            else:
                print('\'%s\' command expects one argument'%args[0])
        elif args[0] == "reverse_complement":
            if len(args) == 2:
                print(reverse_complement(args[1]))
            else:
                print('\'%s\' command expects one argument'%args[0])
        elif args[0] == "calc_nbases":
            if len(args) == 2:
                print('nbases =',calc_nbases(args[1]))
            else:
                print('\'%s\' command expects one argument'%args[0])
        elif args[0] == "is_valid":
            if len(args) == 3:
                print(is_valid(args[1],args[2]))
            else:
                print('\'%s\' command expects two argument'%args[0])
        elif args[0] == "filter_nbases":
            if len(args) == 2:
                print(filter_nbases(args[1]))
            else:
                print('\'%s\' command expects one argument'%args[0])
        elif args[0]=="seq_alignment":
            if (len(args) == 2):
                if (len(opts) == 0):
                    online_alignment(args[1])
                else:
                    for opt, arg in opts:
                        if (opt in ['-o']):
                            online_alignment(args[1], opts[0][1])
                        else:
                            print("wrong option")

            elif (len(args)==3 ):
                x=args[1].rfind("fasta")
                if(x==len(args[1])-5):
                   if(len(opts)==0):
                     seq_alignment_files(args[1],args[2])
                   else:
                     for opt,arg in opts:
                         if(opt in ['-o']):
                           seq_alignment_files(args[1],args[2],opts[0][1])
                         else:
                           print("wrong option")
                else:
                    if (len(opts) == 0):
                        seq_alignment(args[1], args[2])
                    else:
                        for opt, arg in opts:
                            if (opt in ['-o']):
                                seq_alignment(args[1], args[2], opts[0][1])
                            else:
                                print("wrong option")
            else:
                print("No alignment function matches this!!")

        elif args[0] == "convert_to_fasta":
            if len(args) == 2:
                print(convert_to_fasta(args[1]))
            else:
                print('\'%s\' command expects file path'%args[0])
        elif args[0] =="merge_fasta":
            if len(args) >= 3:
                if opts!=[]:
                    merge_fasta(*args[1:], outputPath=opts[0][1])
                    print('merging is done in \'%s\''%opts[0][1])
                else:
                    merge_fasta(*args[1:])
            else:
                print('\'%s\' command expects at least two fasta file'%args[0])
        else:
            print("\'%s\' is not recognized as a Biological Command."%args[0])
if __name__ == "__main__":
    main()
