from modules.file_handling import get_args
from modules.file_handling import read_fasta

from modules.nucleotide_entropy import NucleotideEntropy

import matplotlib.pyplot as plt

from itertools import tee


def pairwise(it):
    it = iter(it)
    while True:
        try:
            yield next(it), next(it)
        except StopIteration:
            return

def downsample(old_list, fold: int):
	new_list = []
	for a, b in pairwise(old_list):	
		new_list.append((a + b) / 2)
		new_list.append((a + b) / 2)
	fold -= 1
	fold = 1
	if fold > 1:
		new_list = downsample(new_list, fold)
	return new_list

#----------------------------command line arguements--------------------------#

args = get_args()

#----------------------------------file input---------------------------------#

contig_dict = read_fasta(args.infile)

#------------------------------------stuff------------------------------------#
fig, ax = plt.subplots()


for id in contig_dict:
	contig_entropy = NucleotideEntropy(contig_dict[id])

	
	data  = contig_entropy.get()
	for frame in [1,2,3]:
		ax.plot(list(data[frame])[0:1000], label='frame ' + str(frame))

a = ax.get_xticks().tolist()
a = [x * 3 for x in a]
ax.set_xticklabels(a)

ax.legend(loc="upper left")
fig.set_size_inches(20, 5)
fig.savefig('figure.png', dpi=300)
plt.show()
