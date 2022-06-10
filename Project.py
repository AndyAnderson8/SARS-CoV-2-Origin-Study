import os #module for file access
import jellyfish #module for mutation algorithm
import matplotlib.pyplot as plt #module for plotting

compareFileName = "SARS-CoV-2_COMPARE.fasta"
fastaFolderDir = "FASTAs"

aminoAcidDict = { #dictonary turning RNA codons to corresponding amino acid
  "UUU": "F", "UUC": "F",
  "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
  "AUU": "I", "AUC": "I", "AUA": "I",
  "AUG": "M",
  "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
  "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
  "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
  "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
  "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
  "UAU": "Y", "UAC": "Y",
  "UAA": "!", "UAG": "!", "UGA": "!", #stop codon
  "CAU": "H", "CAC": "H",
  "CAA": "Q", "CAG": "Q",
  "AAU": "N", "AAC": "N",
  "AAA": "K", "AAG": "K",
  "GAU": "D", "GAC": "D",
  "GAA": "E", "GAG": "E",
  "UGU": "C", "UGC": "C",
  "UGG": "W",
  "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
  "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def loadFastas(directory):
  '''Takes in all fasta files from specified folder and creates array of the sequences contained within them'''
  fastaFiles = os.listdir(directory)
  fastaSequences = [] #holds all fasta sequences
  for file in fastaFiles:
    openFile = open(directory + "/" + file, "r")
    lineArray = openFile.readlines()
    lineArray.pop(0) #remove header line, not used
    sequence = ""
    for line in lineArray:
      line = line.replace("T", "U").upper() #makes sure it's RNA, not DNA
      sequence += line.strip() #removes irregular formatting
    fastaSequences.append([file, sequence])
  return fastaSequences
  
def translate(rnaSequence):
  '''Converts RNA sequence to the amino acid sequence it codes for'''
  aminoAcidSequence = ""
  started = False #checks for start codon before translating
  while len(rnaSequence) >= 3: #checks if there is more codons to translate
    codon = rnaSequence[0:3]
    if "R" in codon or "Y" in codon or "K" in codon or "M" in codon or "S" in codon or "W" in codon or "B" in codon or "D" in codon or "H" in codon or "V" in codon or "N" in codon or "-" in codon:
      aminoAcid = "X"
    else:
      aminoAcid = aminoAcidDict[codon]
      if aminoAcid == "M": #start codon
        started = True
      elif aminoAcid == "!": #stop codon checking is only really important for more than one protein
        started = False
    if started:
      aminoAcidSequence += aminoAcid
      rnaSequence = rnaSequence[3:len(rnaSequence)] #remove codon from rnaSequence
    else:
      rnaSequence = rnaSequence[1:len(rnaSequence)] #remove single base pair, continue with open reading frame
  return aminoAcidSequence

rnaSequences = loadFastas(fastaFolderDir) #relative path to folder with FASTA files
aminoAcidSequences = [] #array of all varients amino acids
baseAminoAcidSequence = [] #stores sequence of alpha varient to compare to

for rnaSequence in rnaSequences: #translate all RNA sequences
  rnaSequence[1] = translate(rnaSequence[1])
  if rnaSequence[0] == compareFileName:
    baseAminoAcidSequence = rnaSequence
  else:
    aminoAcidSequences.append(rnaSequence)

for aminoAcidSequence in aminoAcidSequences: #count mutations in all sequences
  mutationCount = jellyfish.levenshtein_distance(baseAminoAcidSequence[1], aminoAcidSequence[1]) #how many mutations
  print("Mutations between " + baseAminoAcidSequence[0] + " and " + aminoAcidSequence[0] + ": " + str(mutationCount))