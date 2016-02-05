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
#   Samelson_HW2.py -f sample_5mer -l 5
#
##############################################################################
#
#   Important assumptions:
#
#   1- Assuming FASTA file is small enough that a small personal computer can run it
#
#   2- Assuming sequences are legitimate and are not filled with random junk
#
##############################################################################
###############################################################################
#
#   Speed and memory limitations
#
#   O(n^2)
#
#   I don't have proof, however, the algorithm is a nested for loop with a single
#   calculation inside.  This algorithm is not the most efficient, but it is not
#   the worst case either.  There is an algorithm similar to mine called, the Nussinov
#   algorithm that runs O(n^2), so I have sufficient reason to believe my algorithm runs
#   in similar time.
#
##############################################################################
#
#   Outline:
#   I - Parsing command line arguments
#   II - Finding most common Kmer
#   III - Printing the top 5 most common Kmers
#   IV - Reading FASTA file
#
##############################################################################
#
#  References/Collaborators:
#
#      1.  Collaborated with Ali Hakimi and Chris Gray
#
#      2.  Commenting style formatted after Korshoj
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
#   Properly assures user inputs correct values and throws errors if inputs are not correct
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

    filename = None   # Name of FASTA file
    chromo = None     # Chromosome name in the FASTA file
    kmer_size = None  # size of k-mer to search for

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

##############################################################################
#
#   II- Finding most common kmer of each
#
#   1- Function finds all unique kmers of the given size and counts their frequencies in a dictionary
#
#   2. Function does this for every sequence in the file
#
#   3- It then tallies up all the totals and prints the top 5 most common kmers of the given size
#
#   for every sequence in the file:
#        for i to end of each sequence:
#               find unique kmers
#        populate second dictionary with number of times kmer appears
#        reset first dictionary
#   Reverse the final dictionary, most common to least
#   Return final dictionary
##############################################################################

def find_most_common(sequenceArray, size):


    #Create dictionaries to store unique kmers and frequencies
    freq = defaultdict(int)
    finaldict = defaultdict(int)

    # Loop through every sequence in the file
    for s in sequenceArray:
        #print "Unique Kmers" + "\n"
        for i in range(len(s) - size + 1):
            freq[s[i:i + size]] += 1       #find unique kmers and add them to dictionary

        for key in freq:                   #populate second dictionary with number of times kmer appears in file
            if key in finaldict:
                finaldict[key] += 1
            else:
                finaldict[key] = 1

        freq = defaultdict(int)            #reset the first dictionary


    finaldict = finaldict.items()
    finaldict.sort(key = lambda item: item[1], reverse=True)  #reverse sort the final dict, highest to lowest

    return finaldict




##############################################################################
#
#   III- Printing the top 5 most common Kmers
#
#   1- Function takes the dictionary, takes the top 5, and outputs them to screen
#
##############################################################################

def print_top_5(finaldict, size):
    orig_stdout = sys.stdout
    finaldict = finaldict[:5]
    print "################################################"
    print "The top 5 most common " + str(size) +"mers that appear are..."
    print "################################################"
    for i in finaldict:
        print i[0] + " " + str(i[1])

    f = open("SAMPLE.OUT", 'w')                    # output to SAMPLE.OUT
    sys.stdout = f
    print "################################################"
    print "The top 5 most common " + str(size) +"mers that appear are..."
    print "################################################"
    for i in finaldict:
        print i[0] + " " + str(i[1])
    sys.stdout = orig_stdout
    f.close()



##############################################################################
#
#   IV- Reading the FASTA file
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
        if line.startswith('>') == False:                   #Add to sequence if the line doesn't start with '<'
            sequence += line.strip()
            if i == len(lines) -1:
                sequences.append(sequence)
                break
        else:
            if sequence:
                sequences.append(sequence)
            sequence = ""

    return sequences                                        #Return an array of all the sequences in the file



if __name__ == '__main__':
    [filename, kmer_size] = parsing()                  # Parse the command line return filename and kmer size to search
    kmer_size = int (kmer_size)
    seq = read_fasta_file(filename, kmer_size)         # Read  through the file and return an array of sequences
    finaldict = find_most_common(seq, kmer_size)       # Find most common kmers of given size
    print_top_5(finaldict, kmer_size)                  # Output results to command line


