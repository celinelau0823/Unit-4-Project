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

rnaSequence = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
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
  return aminoacids
 
print(rna2protein(rnaSequence, rnacd.RNAcodonDict))

def protMass(aminoacids):
    weightMass = 0 #starting sum = 0
    for aacid in aminoacids:
      weightMass += massTable.mIsoTable[aacid]
    return weightMass

rosalind_10 = 'ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG'
rosalind_12 = 'ATCGGTCGAA'
rosalind_15 = 'ATCGGTCGAGCGTGT'

def mRNA(dnaSeq, intronOne, intronTwo): #RNA SPLICE
  if intronOne in dnaSeq:
    while intronOne in dnaSeq:
      dnaSeq.remove(intronOne)
  
  if intronTwo in dnaSeq:
    while intronTwo in dnaSeq:
      dnaSeq.remove(intronTwo)

  tempRNA = dna2rna(dnaSeq)
  splicedRNA = rna2protein(tempRNA)
  return splicedRNA

print(mRNA(rosalind_10,rosalind_12,rosalind_15))




