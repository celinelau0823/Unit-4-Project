###Celine Lau
###Lab 4: Project Rosalind

import codon as d# this is my file with your codon dictionaries!
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

dnaSeq = 'GATGGAACTTGACTACGTAAATT'

def dna2rna(dnaSeq):
 '''translates DNA into RNA'''
 rnaSeq = ""
 for T in dnaSeq:
     rnaSeq = dnaSeq.replace('T', 'U') 
 print(rnaSeq)
 return rnaSeq


#DNA to protein code from lab 2  
#uses the codon table to translate groups of 3 base pairs into a protein
def dna2protein(dnaSeq, codonTableD):
  '''Given a dnaSeq, output corresponding aa chain.

  Input: str of dnaSeq with len%3=0
         dict, codonTableD with keys=codon and vals=aa
  Output: str of amino acid sequence
  '''
  aminoacids = ""
  if dnaSeq == "":
    return "" 
  elif len(dnaSeq) % 3 != 0: #if there are extra nucleotides that cannot be bunched into 3, it is cut off
    dnaSeq = dnaSeq[:-(len(dnaSeq) %3)]

  #converts 3 amino acids into its respective protein, and then adds it to the protein chain
  for i in range(0, len(dnaSeq), 3):
    aminoacids += codonTableD[dnaSeq[i:i + 3]]
  return aminoacids

#reformatted version of reverse strand maker from lab 2
#reverses order of original strand and creates compliment
def rStrand(dnaSeq):
  while len(dnaSeq) % 3 != 0:
    dnaSeq = dnaSeq[:len(dnaSeq)-1]
  r = dnaSeq[::-1] #indexes dnaSeq backwards

  rComp = r.replace('A', 't').replace('T', 'a').replace('C','g').replace('G','c') #creates opposite strand of DNA
  rComp = rComp.upper() #reformats the dna strand
  return rComp

