import sys, getopt
from collections import Counter
from collections import defaultdict
from operator    import itemgetter

###########################################################################
#  Main function.  Takes in command line options, opens the file, calculates GC content, finds the K-mers
#  and outputs relevant information to "output.gff3"
#
#  Assumption #1:  User uses command prompt correctly using the format test.py -f <filename> -c <chromo name> -k <k-mer>
#  Assumption #2:  User gives a legitimate chromosome with proper nucleotides, not random letters or characters
#  Assumption #3:  User's chromosome name exists in file and K-mer is workable
#########################################################################



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



def most_common_kmer(seq, k_size):
    accumulator = Counter(seq)
    #for length in range (k_size, len(seq) + 1):
        #for start in range (len(seq) - length):
            #accumulator[seq[start:start+length]] += 1
    #print accumulator.most_common(5)

def find_most_common(s, X):
    freq = defaultdict(int)
    for i in range(len(s) - X + 1):
        freq[s[i:i+X]] += 1
    freq = freq.items()
    freq.sort(key = lambda item: item[1], reverse=True)
    top_five = freq[:5]
    print top_five
    print "\n"


def read_fasta_file(filename, kmer_size):
    file_open = open(filename,"r")
    line = file_open.readline()
    if not line.startswith(">"):
        raise TypeError("Not a FASTA file: %r" % line)
    title = line[1:].rstrip()
    seq = []
    while 1:
        line = file_open.readline().rstrip()
        if line == "":
            break
        seq.append(line)
    seq = "".join(seq)
    seq.replace("\n", "")
    print("sequence name is " + title +  "\n\n")
    print ""
    print("The top 5 most common " + str(kmer_size) + "mers  are..." + "\n")
    return seq


if __name__ == '__main__':
    [filename, kmer_size] = parsing()
    kmer_size = int (kmer_size)
    seq = read_fasta_file(filename, kmer_size)
    #most_common_kmer(seq, kmer_size)
    find_most_common(seq, kmer_size)