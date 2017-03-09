#!/usr/bin/python
# Importing libraries
import sys
import numpy as np

# Read Fasta file in a dictionary  
def readFASTA(input_fasta):
  fh = open(input_fasta, 'r')
  seq_dict = {}
  for line in fh:
    line = line.rstrip('\n')
    if '>' in line:
	  name = line
	  seq_dict[name] = ''
    else:
	  seq_dict[name] += line
  return seq_dict

#Build markov_chain to do further computation
def build_markov_chain(sequence_dict):
  markov_chain = {}
  sequence = ''
  for element in sequence_dict.values():
	sequence += element
  #print sequence
  nucl = list(sequence)
  AA, AT, AC, AG, TA, TT, TC, TG, CA, CT, CC, CG, GA, GT, GC, GG = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  for i in range(len(nucl)-1):
	if nucl[i] == 'A' and nucl[i+1] == 'A':
	  AA += 1
	elif nucl[i] == 'A' and nucl[i+1] == 'T':
	  AT += 1
	elif nucl[i] == 'A' and nucl[i+1] == 'C':
	  AC += 1
	elif nucl[i] == 'A' and nucl[i+1] == 'G':
	  AG += 1
	elif nucl[i] == 'T' and nucl[i+1] == 'A':
	  TA += 1
	elif nucl[i] == 'T' and nucl[i+1] == 'T':
	  TT += 1
	elif nucl[i] == 'T' and nucl[i+1] == 'C':
	  TC += 1
	elif nucl[i] == 'T' and nucl[i+1] == 'G':
	  TG += 1
	elif nucl[i] == 'C' and nucl[i+1] == 'A':
	  CA += 1
	elif nucl[i] == 'C' and nucl[i+1] == 'T':
	  CT += 1
	elif nucl[i] == 'C' and nucl[i+1] == 'C':
	  CC += 1
	elif nucl[i] == 'C' and nucl[i+1] == 'G':
	  CG += 1
	elif nucl[i] == 'G' and nucl[i+1] == 'A':
	  GA += 1
	elif nucl[i] == 'G' and nucl[i+1] == 'T':
	  GT += 1
	elif nucl[i] == 'G' and nucl[i+1] == 'C':
	  GC += 1
	elif nucl[i] == 'G' and nucl[i+1] == 'G':
	  GG += 1
  markov_chain['AA'] = float(AA) / (AA+AT+AC+AG)  
  markov_chain['AT'] = float(AT) / (AA+AT+AC+AG)  
  markov_chain['AC'] = float(AC) / (AA+AT+AC+AG)  
  markov_chain['AG'] = float(AG) / (AA+AT+AC+AG)  
  markov_chain['TA'] = float(TA) / (TA+TT+TC+TG)  
  markov_chain['TT'] = float(TT) / (TA+TT+TC+TG)  
  markov_chain['TC'] = float(TC) / (TA+TT+TC+TG)  
  markov_chain['TG'] = float(TG) / (TA+TT+TC+TG)  
  markov_chain['CA'] = float(CA) / (CA+CT+CC+CG)  
  markov_chain['CT'] = float(CT) / (CA+CT+CC+CG)  
  markov_chain['CC'] = float(CC) / (CA+CT+CC+CG)  
  markov_chain['CG'] = float(CG) / (CA+CT+CC+CG)  
  markov_chain['GA'] = float(GA) / (GA+GT+GC+GG)  
  markov_chain['GT'] = float(GT) / (GA+GT+GC+GG)  
  markov_chain['GC'] = float(GC) / (GA+GT+GC+GG)  
  markov_chain['GG'] = float(GG) / (GA+GT+GC+GG) 

  return markov_chain

#Calculate likelihood   
def cal_likelihood(sequence):
  seq = list(sequence)
  background = 2.5 
  markov = 2.5
  for i in xrange(1, len(seq)):
    pair = seq[i-1]+seq[i]
    background *= 1.25
    markov = markov *markov_chain[pair] * 5
  likelihood = np.log10(markov/background)
  #print background
  #print markov
  #print likelihood
  return likelihood

 #Translate DNA sequence to a protein sequence
def translation(sequence):
    codon_table = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
    'TAT':'Y', 'TAC':'Y', 'TAA':'_', 'TAG':'_', 'TGT':'C', 'TGC':'C', 'TGA':'_', 'TGG':'W',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCG':'A', 'GCA':'A',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    }
        
    protein_sequence = ''    
    for n in range(0, len(sequence), 3): #translate the sequence from position0
        if codon_table.has_key(sequence[n:n+3]) == True: #make sure it doesn't return errors if there are less than 3 bases left
            protein_sequence += codon_table[sequence[n:n+3]] #translate based on the dictionary above
            
    return protein_sequence  

#A format for output file
def print_format(input_dict, ouput_gene, output_protein):
  output_dict = {}
  input_count = 0
  for key in input_dict.keys():
    input_count += 1
    likelihood = cal_likelihood(input_dict[key])
    output_dict[key] = likelihood

  fg = open(ouput_gene, 'w')
  fp = open(output_protein, 'w')
  
  print "IS gene and it's likelihood:"
  predicted_count = 0
  for key in output_dict.keys():
	if output_dict[key] > 0:
		predicted_count += 1
		print key
		print output_dict[key]
		fg.write(">" + key + "\n")
		fg.write(input_dict[key] + "\n")
		fp.write(">" + key + "\n")
		protein_seq = translation(input_dict[key])
		fp.write(protein_seq + "\n")
  accuracy = float(predicted_count)/input_count*100
  print "The accuracy of my prediction is %0.4f percent" %accuracy

#Main Function
if __name__ == "__main__":
  if len(sys.argv) != 3:
	print "USAGE: python ECgnfinder_mc.py -i [INPUT_FILE]"
	exit(-1)

  if sys.argv[1] == "-i":
	if sys.argv[2].endswith('.fasta'):
		input_file = sys.argv[2]
	else:
		print "The input file should be in fasta format"
		exit(-1)
  else:
		print "USAGE: python ECgnfinder_mc.py -i [INPUT_FILE]"
		exit(-1)
		
  sequence_dict = readFASTA("gene_1000.fasta")
  markov_chain = build_markov_chain(sequence_dict)
  input_dict = readFASTA(input_file)
  print_format(input_dict, "output_gene.fasta", "output_protein.fasta")
