##############################################################################
##############################################################################
#
#   File     : Samelson_HW2.py
#
#   Purpose  : This program is designed for finding the most common k_mer of given
#              length.  It searches through a given FASTA file and outputs the top
#              5 most common k_mers of length k.
#
#   Developer : Lincoln Samelson
#               CSCI 4314, Spring 2016, HW2
#
##############################################################################
#
#   Sample command line arguments to run program:
#   Samelson_HW2.py -f YeastGenome.fasta -l 7
#
##############################################################################
#
#   Important assumptions:
#
#   1-
#
##############################################################################
#
#   Outline:
#   I - Parsing command line arguments
#   II - Finding most common Kmer
#   III - Reading FASTA file
#
##############################################################################

##############################################################################
#
#   I - Parsing command line arguments into variables
#
#   User must enter arguments as denoted (order does not matter):
#   -f <FASTA filename>
#   -l k_mer size to search for
#
#   The 'getopt' and 'sys' sets are imported and used to parse incoming
#   arguments. A help function is called if the user does not enter arguments
#   correctly.
#
##############################################################################


import sys, getopt
from collections import Counter
from collections import defaultdict
from operator    import itemgetter

def user_help():

    sys.stdout.write("\n\n\nCommand line arguments must be entered in the following format, using all specifiers:\n\n")
    sys.stdout.write("  -f   <FASTA filename>            \n")
    sys.stdout.write("  -l   <size of k-mer to search>  \n\n")


def parsing():

    filename = None # Name of FASTA file
    chromo = None # Chromosome name in the FASTA file
    kmer_size = None # size of k-mer to search for

    # Getting the arguments from the user
    try:
        options, remainder = getopt.getopt(sys.argv[1:],"f:l:", ["help"])

    # An error ends the program and displays the help comments
    except getopt.GetoptError:
        sys.stdout.write("\n\n\nArguments were not entered correctly.")
        user_help()
        sys.exit(1)

    # Processing the list of arguments from 'getopt'
    for op, value in options:
        if op=="-f":
            filename = value
            print "\n"
        elif op=="-l":
            kmer_size = value
            kmer_size = int(kmer_size)
            if kmer_size > 8 or kmer_size < 3:
                print "Sorry!  Kmer size must be between 3 and 8.  Try again."
                sys.exit(1)
            print "\n"
        elif op=="--help":
            user_help()
        else:
            sys.stdout.write("Unhandled argument: [%s][%s]" % (op))

    # The user must enter all arguments or the program is stopped
    if filename == None  or kmer_size == None:
        sys.stdout.write("\n\n\nNot all arguments were specified.")
        user_help()
        sys.exit(1)

    return filename, kmer_size



# def most_common_kmer(seq, k_size):
#     accumulator = Counter(seq)
#     #for length in range (k_size, len(seq) + 1):
#         #for start in range (len(seq) - length):
#             #accumulator[seq[start:start+length]] += 1
#     #print accumulator.most_common(5)

##############################################################################
#
#   II- Finding most common kmer
#
#   1- Program finds all subsequences of given size and counts their frequencies in a dictionary
#
#   2- It then prints the top 5 most common kmer of the given size
##############################################################################

def find_most_common(sequenceArray, size):
    # s = 'hello'
    # sequenceArray = ['hellohellohello', 'skoooooooooooooo']
    #freq = defaultdict(int)
    # currentkmer = {}
    # kmer = ''
    # numChars = 0
    # found = 0
    #
    # for sequence in s:
    #     for char in enumerate(sequence):
    #         index = char[0]
    #         seqLength = len(sequence)
    #         for i in range(index, index + size):
    #             if index <= seqLength - size:
    #                 kmer += sequence[i]
    #
    #         for kmerSequence



    freq = defaultdict(int)
    finaldict = defaultdict(int)
    #print freq
    for s in sequenceArray:
        #print "Unique Kmers" + "\n"
        for i in range(len(s) - size + 1):
            freq[s[i:i + size]] += 1

        A = Counter(freq)
        B = Counter(finaldict)

        finaldict = A + B
        freq = freq.items()
        #print freq
        #print '\n'
        freq = defaultdict(int)
    # print top_five
    # print "\n"
    #print finaldict
    finaldict = finaldict.items()
    finaldict.sort(key = lambda item: item[1], reverse=True)
    #print finaldict
    top_five = finaldict[:5]

    return finaldict

def print_top_5(finaldict):
    finaldict = finaldict[:5]
    print "################################################"
    print "The top 5 most common Kmers in this FASTA File"
    print "################################################"
    for i in finaldict:
        print i[0] + " " + str(i[1])
##############################################################################
#
#   III- Reading the FASTA file
#
#   1- Opens the FASTA File
#
#   2- Gets name of sequence and the actual sequence.
#
#   3- Checks to make sure it is actually a FASTA file
##############################################################################

def read_fasta_file(filename, kmer_size):
    f = open(filename, 'r')
    lines = f.readlines()
    sequence = ''
    sequences = []
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith('>') == False:
            sequence += line.strip()
            if i == len(lines) -1:
                sequences.append(sequence)
                break
        else:
            if sequence:
                sequences.append(sequence)
            sequence = ""
    #print sequences
    #print sequence
    return sequences
    # file_open = open(filename,"r")
    # line = file_open.readline()
    # if not line.startswith(">"):
    #     raise TypeError("Not a FASTA file: %r" % line)
    # title = line[1:].rstrip()
    # seq = []
    # while 1:
    #     line = file_open.readline().rstrip()
    #     if line == "":
    #         break
    #     seq.append(line)
    # seq = "".join(seq)
    # seq.replace("\n", "")
    # print("sequence name is " + title +  "\n\n")
    # print ""
    # print("The top 5 most common " + str(kmer_size) + "mers are..." + "\n")
    # return seq


if __name__ == '__main__':
    [filename, kmer_size] = parsing()
    kmer_size = int (kmer_size)
    seq = read_fasta_file(filename, kmer_size)
    finaldict = find_most_common(seq, kmer_size)
    print_top_5(finaldict)