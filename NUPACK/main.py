import numpy as np
import matplotlib.pyplot as plt
from nupack import *
import os
from utils import *
from shutil import copyfile

comFile = 'commands.UTP.dat' # name of command file
#mmdPath = 'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'

#testSeq ='TAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC' #UTP seq


def numbers2letters(sequences): #Tranforming letters to numbers:
    '''
    Converts numerical values to ATGC-format
    :param sequences: numerical DNA sequences to be converted
    :return: DNA sequences in ATGC format
    '''
    if type(sequences) != np.ndarray:
        sequences = np.asarray(sequences)

    my_seq=["" for x in range(len(sequences))]
    row=0
    for j in range(len(sequences)):
        seq = sequences[j,:]
        assert type(seq) != str, 'Function inputs must be a list of equal length strings'
        for i in range(len(sequences[0])):
            na = seq[i]
            if na==0:
                my_seq[row]+='A'
            elif na==1:
                my_seq[row]+='T'
            elif na==2:
                my_seq[row]+='C'
            elif na==3:
                my_seq[row]+='G'
        row+=1
    return my_seq

'''
fastaFile = 'sequencesUTPtruncated24to85.txt'
f = open(fastaFile,'r')
seqs = f.read()
f.close()
seqs = seqs.split('\n')[0]
seqs = seqs.split(',')
'''
seq = numbers2letters(np.random.randint(0,4,size=(10,100)))[0]
seqs = [seq,seq]
#seqs = ['TATGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC','TATGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC']#['TAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC','TTAGATCCGCATGAGGCTCGATCTGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC']



def getSecondaryStructure(testSeq):
    #testSeq = 'GCGTTCATTCCGC' # nice hairpin
    #testSeq = 'GTGCCTCTTGAGGGTAATGGTTTGATTAACCGGAAGTGAT' # random
    A = Strand(testSeq,name='A') # good UTP binding strand
    #A = Strand(makeFasta(40),name='A')
    #B = Strand(makeFasta(40),name='B')
    #C = Strand(makeFasta(40),name='C')

    #comp = Complex([A,B,C], name='ABC')
    #set1 = ComplexSet(strands=[A,B,C], complexes=SetSpec(max_size=1, include=[comp]))
    comp = Complex([A], name='AA')
    set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1, include=[comp]))


    model1 = Model(material='dna', celsius = 37.0, sodium=0.05)
    results = complex_analysis(set1, model=model1, compute=['pfunc', 'pairs', 'subopt', 'mfe'])
    cout = results[comp]

    print('\nMFE proxy structures(s) for set')
    for i, s, in enumerate(cout.mfe):
        print('    %2d: %s (%.2f kcal/mol)' % (i, s.structure, s.energy))

    #print('\nSuboptimal proxy structures for set')
    #for i, s in enumerate(cout.subopt):
    #    print('    %2d: %s (%.2f kcal/mol)' % (i, s.structure, s.energy))


    # extract list of most probable base pairs
    pairMat = cout.pairs.to_array()
    pairList = []
    paired = []
    for i in range(len(pairMat)): # this whole thing doesn't seem to work properly to me
        bestMatch = np.argmax(pairMat[i,:])
        if True:#pairMat[i,bestMatch] > 0.5: # if we're relatively confident about this secondary structure feature
            if bestMatch != i: # if the best match is itself, it means it does not bind
                if (not bestMatch in paired) and (not i in paired):# check for duplicates
                    pairList.append([i + 1,bestMatch + 1])
                    paired.append(i)
                    paired.append(bestMatch)

    return s.structure, pairList


# get secondary structure
pairLists = []
ssStrings = []
for i in range(len(seqs)):
    string, pairs = getSecondaryStructure(seqs[i])
    ssStrings.append(' %s' % string)
    ssStrings[-1]=ssStrings[-1][1:]
    pairLists.append(pairs) # something wrong with these

outdir = {}
outdir['sequences'] = seqs
outdir['strings'] = ssStrings
outdir['pairs'] = pairLists


def ssToList(ssString):
    pairList = []
    paired = []
    for i in range(len(ssString)):
        if ssString[i] == '(': # if it's paired
            counter = 0
            for j in range(1,len(ssString[i:])): # look for the thing paired to it
                if ssString[i+j] == '(':
                    counter += 1
                if ssString[i+j] == ')':
                    if counter > 0:
                        counter -= 1
                    elif counter == 0:
                        # check for duplicates
                        if (not i in paired) and (not i+j in paired):  # check for duplicates
                            paired.append(i)
                            paired.append(i+j)
                            pairList.append([i + 1,i+j + 1]) #make pair list in 1-n basis
                            break

    return pairList

'''MMB
# write pair list as forces to the MMB command file
copyfile('commands.template.dat', comFile) # make command file
replaceText(comFile,'SEQUENCE',testSeq)

baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'
lineNum = findLine(comFile, baseString) # find the line number to start enumerating base pairs

for pair in pairList:
    filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(pair[0],pair[1])
    addLine(comFile, filledString, lineNum + 1)

#os.system(mmdPath + ' -c ' + comFile)
'''

''' doesn't work in wsl
plt.imshow(cout.pairs.to_array())
plt.xlabel('Base index')
plt.ylabel('Base index')
plt.title('Pair probabilities for set')
plt.colorbar()
plt.clim(0, 1)
'''