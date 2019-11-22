from modules.file_handling import get_args
from modules.file_handling import read_fasta

from modules.nucleotide_entropy import NucleotideEntropy

import matplotlib.pyplot as plt



#----------------------------command line arguements--------------------------#

args = get_args()

#----------------------------------file input---------------------------------#

contig_dict = read_fasta(args.infile)

#------------------------------------stuff------------------------------------#
fig, ax = plt.subplots()


for id in contig_dict:
	contig_entropy = NucleotideEntropy(contig_dict[id])
	ax.plot(list(contig_entropy.entropy_at()[1])[50:2000], label='frame 1') 
	ax.plot(list(contig_entropy.entropy_at()[2])[50:2000], label='frame 2') 
	ax.plot(list(contig_entropy.entropy_at()[3])[50:2000], label='frame 3') 
ax.legend(loc="upper left")
fig.set_size_inches(20, 5)
fig.savefig('test.png', dpi=300)
plt.show()
