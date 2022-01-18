###Celine Lau
###Lab 4: Project Rosalind

import rnacodon as rnacd #file with rna to codon dicitonary
import mim as massTable 
import numpy as np
import os

#DNA reformat code from lab 2
#allows the dna sequence to be read and interpreted by the other functions
def read_FASTA(file_path): 
  f = open(file_path, 'r')
  flines = f.readlines()
  lineNum = 0
  d = {} #d is the codon dictionary
  for line in flines:
    lineNum += 1
    if line[0] == ">":
      key = line.strip(">").strip("\n")
      d[key] = ""
    else: #update d[key] so its value includes my current, non-key, sequence line
      line = line.strip("\n")
      d[key] = d[key]+line
  return(d)

def dna2rna(dnaSeq):
 '''translates DNA into RNA
 input: DNA sequence string
 output: RNA sequence string'''
 rnaSeq = ""
 rnaSeq = dnaSeq.replace('T', 'U')  
 return rnaSeq


def rStrand(dnaSeq):
  '''reformatted version of reverse strand maker from lab 2
 reverses order of original strand and creates compliment
 input: DNA sequence string
 output: complementary DNA sequence string'''
  while len(dnaSeq) % 3 != 0:
    dnaSeq = dnaSeq[:len(dnaSeq)-1]
  r = dnaSeq[::-1] #indexes dnaSeq backwards

  rComp = r.replace('A', 't').replace('T', 'a').replace('C','g').replace('G','c') #creates opposite strand of DNA
  rComp = rComp.upper() #reformats the dna strand
  return rComp


#DNA to protein code from lab 2  
#uses the codon table to translate groups of 3 base pairs into a protein
def rna2protein(rnaSeq, codonTableD):
  '''Given a dnaSeq, output corresponding aa chain.

  Input: str of dnaSeq with len%3=0
         dict, codonTableD with keys=codon and vals=aa
  Output: str of amino acid sequence
  '''
  aminoacids = ""
  if rnaSeq == "":
    return "" 
  elif len(rnaSeq) % 3 != 0: #if there are extra nucleotides that cannot be bunched into 3, it is cut off
    rnaSeq = rnaSeq[:-(len(rnaSeq) %3)]

  #converts 3 amino acids into its respective protein, and then adds it to the protein chain
  for i in range(0, len(rnaSeq), 3):
    aminoacids += codonTableD[rnaSeq[i:i + 3]]
  aminoacids = aminoacids.strip('*')
  return aminoacids
 

def protMass(aminoacids):
 '''finds the protein mass from a weighted mass index.
 consult monoisotopic mass index for precise numbers
 input: protein string (string of amino acids)
 output: the weighted mass of all the amino acids summed.'''
 weightMass = 0 #starting sum = 0
 for aacid in aminoacids:
    weightMass += massTable.mIsoTable[aacid] 
 return weightMass


def mRNA(dnaSeq, intronList): 
  '''dna splicing: acts as splicosomes and removes introns from a 
  given DNA sequence. 
  input: a dna sequence and a list of introns
  output: a string of extrons'''
  for intron in intronList:
    if intron in dnaSeq: #only starts if the intron is in the DNA sequence
      while intron in dnaSeq:
        dnaSeq = dnaSeq.replace(intron, '') #replaces intron with a blank, removing it from the string
  
  tempRNA = dna2rna(dnaSeq) #changes spliced sequence into RNA
  extrons = rna2protein(tempRNA, rnacd.RNAcodonDict) #changes new RNA into protein
  extrons = extrons.strip('*') #formatting
  return extrons


def countMutation(dnaSeq, mutatedSeq):
  '''counts point mutation between two DNA Sequences (compares the two)
  input: A normal DNA Sequence and the mutated version of the same sequence
  output: the amount of base pairs that have been mutated'''
  tempMutList = []
  for i in range(len(dnaSeq)):
    if dnaSeq[i] == mutatedSeq[i]:
      tempMutList.append(False) #appends false if not mutated
    elif dnaSeq[i]!= mutatedSeq[i]:
      tempMutList.append(True) #appends true if mutated
  return sum(tempMutList) #gives amount of trues, same as amount of mutations 


def transtranverse(dnaSeq, mutatedSeq):
  '''finds the transition to transversion ratio in point mutations. Refer 
  to explanation in previous docstring for point mutations.

  transition = pyramidine -> pyramidine OR purine -> purine (ie. C->T, A->G, and vice versa of those two)
  transversion = pyramidine -> purine OR vice versa

  input: A normal DNA Sequence and the mutated version of the same sequence
  output: transition to transversion ratio'''
  tempDNAList = []
  for i in range(len(dnaSeq)):
    if dnaSeq[i] != mutatedSeq[i]: #only enter if the DNA has been mutated
      if dnaSeq[i] == 'G' or dnaSeq[i] == 'A':
        if mutatedSeq[i] == 'G' or mutatedSeq[i] == 'A':
          tempDNAList.append(True) #appends true if transition
        else:
          tempDNAList.append(False) #appends false if transversion
      if dnaSeq[i] == 'C' or dnaSeq[i] == 'T':
        if mutatedSeq[i] == 'C' or mutatedSeq[i] == 'T':
          tempDNAList.append(True)
        else:
          tempDNAList.append(False)
  transition = tempDNAList.count(True) #sums up amount of trues in list
  transversion = len(tempDNAList) - transition #total mutations - # of transitions
  return transition/transversion

