from modules.file_handling import get_args
from modules.file_handling import read_fasta
from modules.file_handling import read_gff

from modules.nucleotide_entropy import NucleotideEntropy

import matplotlib.pyplot as plt

from itertools import tee

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

def minimum_frame(l):
	lowest = float("inf")
	frame = None
	for i, f in enumerate([1, -1, 2, -2, 3, -3]):
		if l[i] < lowest:
			lowest = l[i]
			frame = f
	return frame
		

#----------------------------command line arguements--------------------------#

args = get_args()

#----------------------------------file input---------------------------------#

contig_dict = read_fasta(args.infile)

actual_frame = read_gff(args.infile.replace('fna', 'gff'))

#------------------------------------finding error------------------------------------#

for id in contig_dict:
	contig_entropy = NucleotideEntropy(contig_dict[id])
	print(contig_entropy)
	exit()
	correct = 0
	wrong = 0
	for i, base in enumerate(contig_dict[id], 1):
		position = (i-1)//3
		try:
			#print(i, base, actual_frame.get(i, "NC"), minimum_frame(contig_entropy[position]), contig_entropy[position])
			if actual_frame[i] ==  minimum_frame(contig_entropy[position]):
				correct += 1
			else:
				wrong += 1
		except:
			pass #print("0")
	print(correct / (correct + wrong))

exit()


#------------------------------------plotting------------------------------------#
fig, ax = plt.subplots()


for id in contig_dict:
	contig_entropy = NucleotideEntropy(contig_dict[id])
	
	data, data_rev  = contig_entropy.get()

	#forward
	#ave = [sum(a)/3 for a in zip(data[1], data[2], data[3])] 
	#ax.plot(ave, label='forward average')
	for frame in [1,2,3]:
		ax.plot(list(data[frame]), label='frame ' + str(frame))

	#reverse
	#ave = [sum(a)/3 for a in zip(data_rev[1], data_rev[2], data_rev[3])] 
	#ax.plot(ave, label='reverse average')
	for frame in [1,2,3]:
		ax.plot(list(data_rev[frame]), label='frame -' + str(frame))
plt.xlim([0, 1000])
a = ax.get_xticks().tolist()
ax.xaxis.set_major_locator(plt.MaxNLocator(len(a) * 4))
a = ax.get_xticks().tolist()
a = [x * 3 for x in a]
ax.set_xticklabels(a)

ax.legend(loc="upper left")
fig.set_size_inches(20, 5)
fig.savefig('figure.png', dpi=300)
plt.show()
