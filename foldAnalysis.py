'''
stand-alone analysis script
'''
import numpy as np
import os
from seqfold import dg, fold
from tqdm import tqdm
from sty import fg, bg, ef, rs

fastaFile = 'C:/Users\mikem\OneDrive\McGill_Simine\Aptamers\Data/UTP/sequencesUTPtruncated24to85.txt'
f = open(fastaFile,'r')
seqs = f.read()
f.close()
seqs = seqs.split('\n')[0]
fseqs = seqs.join(',')
seqs = seqs.split(',')

#seqs = ['TAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCACT','TAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCACT'] # known good seq - it's #48 in the seq list

goodSeqs = []
goodSeqs.append('TAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC')
#goodSeqs.append('GAAGAGGCTCGATCTAAGGATGTAATTATTCTTTTCGGTTAGTTATCTTCCTAGTTGACTAGTACATGACCACTTGAAAGGGCGAA')
goodSeqs.append('TACATGAGGATGAGGCTCGATCGGTATCTTTGGCATATCCTATCGTTCTATGTAACAGCCTG')
goodSeqs.append('TATTGTCGTTTGAGGCTCGATCCCGCGGTGGGGGCAGGCATATAGTGCCTAGTCTTGGAGTC')
goodSeqs.append('GTAAATTAGTTGAGGCTCGATCTATACTCTTTGTATCTCGTGATTTCCTGGTCGGATTAGCA')

badSeqs = []
badSeqs.append('ACATGGTTCGTGAGGCTCGATCTAGCGTGTAGTGCGGATATATTTGGATAATTGGCAATCAC')
badSeqs.append('ATTGTGTAACTGAGGCTCGATCACGGTGTCCGCCTTATGCCGGTAAGCAGGTCACGGGTGTC')
#badSeqs.append('CCGGTCCGGCTGGCTCGATCGCGTGAGGTTCACTGGGGGGGAGATACACTTGTGGTTTTATTGACTAGTACATGACCACTTGAAAG')
badSeqs.append('ATGTAGTGTGTGAGGCTCGATCGGATTTCGCAGGGGATGATGTTCCCGTTCTCCAGAGCTGG')
badSeqs.append('ACGGCGCCGCTGAGGCTCGATCTAGCCTTCATCTATGTACTATGGTTAGCGGGTCGATTAGA')
badSeqs.append('GCACGTTCATTGAGGCTCGATCGTCAGTGGGTCGTCCCTCGAATCACATTGGGCACGGCTAA')

aa = 1

# take out pair information from the potts model
pottsList = []
offset = 23
def addPottsItem(list, ind1, ind2, offset):
    list.append([ind1-offset, ind2-offset])
    return list


pottsList = addPottsItem(pottsList, 48, 24, offset)
pottsList = addPottsItem(pottsList, 57, 27, offset)
pottsList = addPottsItem(pottsList, 58, 26, offset)
pottsList = addPottsItem(pottsList, 59, 32, offset)
pottsList = addPottsItem(pottsList, 62, 29, offset)
pottsList = addPottsItem(pottsList, 64, 24, offset)
pottsList = addPottsItem(pottsList, 67, 30, offset)
pottsList = addPottsItem(pottsList, 70, 30, offset)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def countBaseFrequency(sequences):
    return sequences.count("A"), sequences.count("T"), sequences.count("G"), sequences.count("C")


def getSecondaryStructure(sequence):  # used in letters2numbers
    '''
    get the secondary structure for a given sequence
    using seqfold here - identical features are available using nupack, though results are sometimes different
    :param sequence:
    :return: a dot-bracket string and list of paired bases (assuming single-strand DNA aptamer)
    '''
    temperature = 37.0  # celcius
    dg(sequence, temp=temperature)  # get energy of the structure
    # print(round(sum(s.e for s in structs), 2)) # predicted energy of the final structure

    structs = fold(sequence)  # identify structural features
    desc = ["."] * len(sequence)
    pairList = []
    for s in structs:
        pairList.append(s.ij[0])
        pairList[-1]  # list of bound pairs indexed from 1
        if len(s.ij) == 1:
            i, j = s.ij[0]
            desc[i] = "("
            desc[j] = ")"

    ssString = "".join(desc)
    pairList = np.asarray(pairList) + 1

    return ssString, pairList


def letters2numbers(sequences): #Tranforming letters to numbers:
    '''
    Converts ATCG sequences to numerical values
    :param sequences: ATCG-format DNA sequences to be converted
    :return: DNA sequences in 1234 format
    '''

    assert type(sequences) == list, 'Function inputs must be a list'

    my_seq=np.zeros((len(sequences), len(sequences[0])))
    row=0

    for seq in sequences:
        assert (type(seq) == str) and (len(seq) == my_seq.shape[1]), 'Function inputs must be a list of equal length strings'
        col=0
        for na in seq:
            if (na=="a") or (na == "A"):
                my_seq[row,col]=0
            elif (na=="u") or (na == "U") or (na=="t") or (na == "T"):
                my_seq[row,col]=1
            elif (na=="c") or (na == "C"):
                my_seq[row,col]=2
            elif (na=="g") or (na == "G"):
                my_seq[row,col]=3
            col+=1
        row+=1

    return my_seq

# get secondary structure
pairLists = []
ssStrings = []
for i in tqdm(range(len(seqs))):
    string, pairs = getSecondaryStructure(seqs[i])
    ssStrings.append(string)
    pairLists.append(pairs)

stringArray = np.asarray(ssStrings)
seqsArray = np.asarray(seqs)
sortedSS = np.sort(stringArray)
sortedSeqs = seqsArray[np.argsort(stringArray)]

#''' print sorted sequences with overlaid secondary structure
#for i in range(len(seqs)):
#    print(i)
#    print('SSString: ' + sortedSS[i].replace('.','-')[0:9] + fg.red + sortedSS[i].replace('.','-')[9:25] + fg.rs + sortedSS[i].replace('.','-')[25:48] + fg.blue + sortedSS[i].replace('.','-')[48:] + fg.rs)
#    print('Sequence: ' + sortedSeqs[i].replace('.','-')[0:9] + fg.red + sortedSeqs[i].replace('.','-')[9:25] + fg.rs + sortedSeqs[i].replace('.','-')[25:48] + fg.blue + sortedSeqs[i].replace('.','-')[48:] + fg.rs)
#'''

# load up nupack's opinions
nuOuts = np.load('C:/Users\mikem\OneDrive\McGill_Simine\Aptamers\Data/UTP/UTP_fold_nupack.npy',allow_pickle=True)
nuOuts = nuOuts.item()
nuStrings = nuOuts['strings']
nuSeqs = nuOuts['sequences']
nuPairs = nuOuts['pairs']

for i in range(len(seqs)):
    print(i)
    #print('NUString: ' + nuStrings[i].replace('.','-')[0:9] + fg.red + nuStrings[i].replace('.','-')[9:25] + fg.rs + nuStrings[i].replace('.','-')[25:48] + fg.blue + nuStrings[i].replace('.','-')[48:] + fg.rs)
    print('SSString: ' + ssStrings[i].replace('.','-')[0:9] + fg.red + ssStrings[i].replace('.','-')[9:25] + fg.rs + ssStrings[i].replace('.','-')[25:48] + fg.blue + ssStrings[i].replace('.','-')[48:] + fg.rs)
    print('Sequence: ' + seqs[i].replace('.','-')[0:9] + fg.red + seqs[i].replace('.','-')[9:25] + fg.rs + seqs[i].replace('.','-')[25:48] + fg.blue + seqs[i].replace('.','-')[48:] + fg.rs)    