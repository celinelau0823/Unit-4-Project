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
 '''translates DNA into RNA'''
 rnaSeq = ""
 rnaSeq = dnaSeq.replace('T', 'U')  
 return rnaSeq

#reformatted version of reverse strand maker from lab 2
#reverses order of original strand and creates compliment
def rStrand(dnaSeq):
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
    weightMass = 0 #starting sum = 0
    for aacid in aminoacids:
      weightMass += massTable.mIsoTable[aacid]
    return weightMass


def mRNA(dnaSeq, intronList): #RNA SPLICE
  for intron in intronList:
    if intron in dnaSeq:
      while intron in dnaSeq:
        dnaSeq = dnaSeq.replace(intron, '')
  
  tempRNA = dna2rna(dnaSeq)
  splicedRNA = rna2protein(tempRNA, rnacd.RNAcodonDict)
  splicedRNA = splicedRNA.strip('*')
  return splicedRNA


def countMutation(dnaSeq, mutatedSeq):
  tempMutList = []
  for i in range(len(dnaSeq)):
    if dnaSeq[i] == mutatedSeq[i]:
      tempMutList.append(False)
    elif dnaSeq[i]!= mutatedSeq[i]:
      tempMutList.append(True)
  return sum(tempMutList)


Rosalind_0209 = 'GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT'
Rosalind_2200 = 'TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT'


def transtranverse(dnaSeq, mutatedSeq):
  tempDNAList = []
  for i in range(len(dnaSeq)):
    if dnaSeq[i] != mutatedSeq[i]:
      if dnaSeq[i] == 'G' or dnaSeq[i] == 'A':
        if mutatedSeq[i] == 'G' or mutatedSeq[i] == 'A':
          tempDNAList.append(True)
        else:
          tempDNAList.append(False)
      if dnaSeq[i] == 'C' or dnaSeq[i] == 'T':
        if mutatedSeq[i] == 'C' or mutatedSeq[i] == 'T':
          tempDNAList.append(True)
        else:
          tempDNAList.append(False)
  transition = tempDNAList.count(True)
  transversion = len(tempDNAList) - transition
  return transition/transversion

print(transtranverse(Rosalind_0209,Rosalind_2200))