import os

from modules.file_handling import get_args
from modules.file_handling import read_fasta
from modules.file_handling import read_gff

from modules.nucleotide_entropy import NucleotideEntropy

import matplotlib.pyplot as plt


def min_dir(a, b):
	pass

def min_idx(a, b, c):
	if(a > b):
		if(b > c):
			return 3;
		else:
			return 2;
	else:
		if(a > c):
			return 3;
		else:
			return 1;

def minimum_frame(position, contig_entropy):
	lowest = float("inf")
	frame = None
	for f in [1, -1, 2, -2, 3, -3]:
		if contig_entropy[f][position] < lowest:
			lowest = contig_entropy[f][position]
			frame = f
	return frame
		

#----------------------------command line arguements--------------------------#

args = get_args()

#----------------------------------file input---------------------------------#
contig_dict = read_fasta(args.infile)

actual_frame = read_gff(args.infile.replace('fna', 'gff'))


for id in contig_dict:
	contig_entropy = NucleotideEntropy(contig_dict[id])
#------------------------------------finding error------------------------------------#
	print('FILENAME', 'AT_SKEW', 'TYPE', sep='\t', end='\t')
	for aa in 'ARNDBCEQZGHILKMFPSTWYV':
		print(aa, end='\t')
	print()
	for i, base in enumerate(contig_dict[id][:-2]):
		dna = contig_dict[id].upper()
		atskew = round((dna.count('A') + dna.count('T')) / len(dna), 2)
		#position = (i-1)//3
		frame = (i % 3 ) + 1
		print(os.path.basename(args.infile), atskew, sep='\t', end='\t')
		try:
			if frame == actual_frame[i+1]:
				print('c', end='\t')
			else:
				print('n', end='\t')
		except:
			print('intergenic', end='\t')
		aminoacid_dictionary = contig_entropy.translate_dict(contig_entropy[0][i])
		for aa in 'ARNDBCEQZGHILKMFPSTWYV':
			print(aminoacid_dictionary.get(aa, 0), end='\t')
		print()
		frame = -frame
		print(os.path.basename(args.infile), atskew, sep='\t', end='\t')
		try:
			if frame == actual_frame[i+1]:
				print('c', end='\t')
			else:
				print('n', end='\t')
		except:
			print('i', end='\t')
		aminoacid_dictionary = contig_entropy.translate_dict(contig_entropy.reverse_frequencies(contig_entropy[0][i]))
		for aa in 'ARNDBCEQZGHILKMFPSTWYV':
			print(aminoacid_dictionary.get(aa, 0), end='\t')
		print()
		#print(i, base, actual_frame.get(i, "NC"), minimum_frame(position, contig_entropy))
		#print(position, '+'+str(frame), sep='\t', end='\t')
		#reverse
		#print(position, '-'+str(frame), sep='\t', end='\t')
		#aminoacid_dictionary = contig_entropy.translate_dict(contig_entropy.reverse_frequencies(contig_entropy[0][i]))
		#for aa in 'ARNDBCEQZGHILKMFPSTWYV':
		#	print(aminoacid_dictionary.get(aa, 0), end='\t')
		#try:
		#	if -frame == actual_frame[i]:
		#		print('coding', end='')
		#	else:
		#		print('noncoding', end='')
				
		#except:
		#	print('intergenic', end='')
			#pass #print(i, position, "bad")
		#print()
	exit()
	dna = contig_dict[id].upper()
	atskew = (dna.count('A') + dna.count('T')) / len(dna)
	forward_skew = (list(actual_frame.values()).count(1) + list(actual_frame.values()).count(2) + list(actual_frame.values()).count(3)) / len(actual_frame.values())
	print("AT_SKEW:", atskew, "PERCENT_CORRECT:", correct / (correct + wrong))
	#exit()
#------------------------------------plotting------------------------------------#
	#continue
	fig, ax = plt.subplots(2)
	#forward
	#ave = [sum(a)/3 for a in zip(data[1], data[2], data[3])] 
	#ax.plot(ave, label='forward average')
	for frame in [1,2,3,-1,-2,-3]:
		ax[0].plot(contig_entropy[frame], label='frame ' + str(frame))
	ax[0].set_xlim(0, 3000)
	a = ax[0].get_xticks().tolist()
	a = [x * 3 for x in a]
	ax[0].set_xticklabels(a)
	ax[0].legend(loc="lower left")
	ax[0].title.set_text('Entropy of Frames')
	ax[0].set_ylabel('Entropy')

	ax[1].scatter(list(actual_frame.keys()), list(actual_frame.values()))
	ax[1].set_yticks([1,2,3,-1,-2,-3])
	ax[1].set_xlim(0, 9000)
	#ax[1].title.set_text('Actual Coding Frame')
	ax[1].set_ylabel('Coding Frame')
	ax[1].set_xlabel('Position Along Genome')

	#plt.xlim([0, 1000])
	
#fig.set_size_inches(20, 5)
#fig.savefig('figure.png', dpi=300)
#plt.show()
