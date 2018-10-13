# -*- coding: utf-8 -*-
"""
Nadya Avirianta S
16/394096/PA/17187

Take Home UTS Bioinformatika
HMM for MSA
"""

import numpy as np
#Biopython untuk meload sekuen yang akan dialign
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq

#Biopython untuk membuat HMM
from Bio.HMM import MarkovModel
from Bio.HMM import DynamicProgramming
from Bio.HMM import Trainer
from Bio.HMM import Utilities

#Definisikan state pada HMM
class DNAAlphabet(Alphabet.Alphabet):
    letters = ['A','T','G','C','']

class StateAlphabet(Alphabet.Alphabet):
    letters = ['I', 'J','K','L','D','E','F','M','N','O']
    
#membangun HMM
hmmbuild = MarkovModel.MarkovModelBuilder(StateAlphabet(),DNAAlphabet())

#transisi dari match state
hmmbuild.allow_transition('M','N')
hmmbuild.allow_transition('N','O')
hmmbuild.allow_transition('M','I')
hmmbuild.allow_transition('N','J')
hmmbuild.allow_transition('O','K')
hmmbuild.allow_transition('M','E')
hmmbuild.allow_transition('N','F')

#transisi dari insert state
hmmbuild.allow_transition('L','M')
hmmbuild.allow_transition('I','N')
hmmbuild.allow_transition('J','O')
hmmbuild.allow_transition('I','I')
hmmbuild.allow_transition('J','J')
hmmbuild.allow_transition('K','K')
hmmbuild.allow_transition('L','L')

#transisi dari delete state 
hmmbuild.allow_transition('D','N')
hmmbuild.allow_transition('E','O')
hmmbuild.allow_transition('D','E')
hmmbuild.allow_transition('E','F')

hmmbuild.allow_transition('N','M')
hmmbuild.allow_transition('O','N')

hmmbuild.allow_transition('I','M')
hmmbuild.allow_transition('J','N')
hmmbuild.allow_transition('K','O')

hmmbuild.allow_transition('M','L')
hmmbuild.allow_transition('N','I')
hmmbuild.allow_transition('O','J')


hmmbuild.allow_transition('E','M')
hmmbuild.allow_transition('F','N')
hmmbuild.allow_transition('N','D')
hmmbuild.allow_transition('O','E')
hmmbuild.allow_transition('E','D')
hmmbuild.allow_transition('F','E')

#initial probability
hmmbuild.set_initial_probabilities({'M': 0.8, 'D': 0.2,'L':0})

#probabilitas untuk transisi dari satu state ke state lain
hmmbuild.set_transition_score('M', 'N', 1)
hmmbuild.set_transition_score('N', 'O', .4)
hmmbuild.set_transition_score('N', 'J', .6)
hmmbuild.set_transition_score('J', 'J', (4/6))
hmmbuild.set_transition_score('J', 'O', (2/6))
hmmbuild.set_transition_score('D', 'N', 1)

#emission probability dari state
hmmbuild.set_emission_score('M', 'A', 1)
hmmbuild.set_emission_score('N', 'G', 1)
hmmbuild.set_emission_score('O', 'C', .8)
hmmbuild.set_emission_score('O', 'G', .2)
hmmbuild.set_emission_score('J', 'A', (6/7))
hmmbuild.set_emission_score('J', 'G', (1/7))
hmmbuild.set_emission_score('D', 'A', 0)
hmmbuild.set_emission_score('D', 'T', 0)
hmmbuild.set_emission_score('D', 'G', 0)
hmmbuild.set_emission_score('D', 'C', 0)
hmmbuild.set_emission_score('E', 'A', 0)
hmmbuild.set_emission_score('E', 'T', 0)
hmmbuild.set_emission_score('E', 'G', 0)
hmmbuild.set_emission_score('E', 'C', 0)
hmmbuild.set_emission_score('F', 'A', 0)
hmmbuild.set_emission_score('F', 'T', 0)
hmmbuild.set_emission_score('F', 'G', 0)
hmmbuild.set_emission_score('F', 'C', 0)

#membangun hmm 
hmm = hmmbuild.get_markov_model()

#contoh sequen yang akan dialign
bat = Seq('AGC',generic_dna)
rat = Seq('AGAGC',generic_dna)
cat = Seq('AGAAG',generic_dna)
gnat = Seq('GAAAC',generic_dna)
goat = Seq('AGC',generic_dna)

#list state untuk masing2 sekuen
catst = MutableSeq('MNJJO',StateAlphabet())
batst = MutableSeq('MNO',StateAlphabet())
ratst = MutableSeq('MNJJO',StateAlphabet())
gnatst = MutableSeq('DNJJO',StateAlphabet())
goatst = MutableSeq('MNO',StateAlphabet())

seq = [cat,bat,rat,gnat,goat]
states = [catst,batst,ratst,gnatst,goatst]

#untuk training HMM
trainer = Trainer.KnownStateTrainer(hmm)
for i in range(len(seq)):
    trainingseq = Trainer.TrainingSequence(seq[i],states[i])
    trainedhmm = trainer.train([trainingseq])
    
#mengetest HMM yang sudah di train
teststates = MutableSeq('', StateAlphabet())

for testseq in seq:
    #digunakan Viterbi algorithm
    predictedstates, prob = trainedhmm.viterbi(testseq, StateAlphabet())
    print("Prediction probability: %f" % prob)
    #untuk mengeprint test sequence dan predictednya
    Utilities.pretty_print_prediction(testseq, teststates, predictedstates)
    
#contoh untuk print alignment antara 2 DNA sekuen
testseq='AGAGC'
profile='ATC'
predictedstates, prob = trainedhmm.viterbi(testseq, StateAlphabet())

print('Predicted state:' +predictedstates)
for i in range(len(predictedstates)):
    if (predictedstates[i]=='D') or(predictedstates[i]=='E')or(predictedstates[i]=='F') :
        testseq = testseq[:i] +'-'+testseq[i:]
        
    elif (predictedstates[i] == 'I') or (predictedstates[i] == 'J') or(predictedstates[i] == 'K')or(predictedstates[i] == 'L'):
        profile=profile[:i]+'-'+profile[i:]
    else:
        testseq=testseq

print('MSA:')
print(testseq)
print(profile)
