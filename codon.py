codonDict = {
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'TTA': 'L', 'TTG': 'L',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TTT': 'F', 'TTC': 'F', 'ATG': 'M', 'TGT': 'C', 'TGC': 'C',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'CCT': 'P',
    'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGG': 'W', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'CAT': 'H', 'CAC': 'H', 'GAA': 'E', 'GAG': 'E', 'GAT': 'D', 'GAC': 'D', 'AAA': 'K',
    'AAG': 'K', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'TAA': '*', 'TAG': '*',
    'TGA': '*'
}

aaDict = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 
    'F': ['TTT', 'TTC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 
    'K': ['AAA', 'AAG'], 'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'], 'M': ['ATG'], 'N': ['AAT', 'AAC'], 
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], '*': ['TAA', 'TAG', 'TGA'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC']
}

cDict = {}
for key in codonDict:
    newKey = key.replace("T", "U")
    newVal = codonDict[key]
    cDict[newKey] = newVal
 

