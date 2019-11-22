from modules.file_handling import get_args
from modules.file_handling import read_fasta

from modules.nucleotide_entropy import NucleotideEntropy

from math import copysign



#----------------------------command line arguements--------------------------#

args = get_args()

#----------------------------------file input---------------------------------#

contig_dict = read_fasta(args.infile)

#------------------------------------stuff------------------------------------#

for id in contig_dict:
	contig_entropy = NucleotideEntropy(contig_dict[id])
	for item in zip(contig_entropy.entropy_at()[1], contig_entropy.entropy_at()[2], contig_entropy.entropy_at()[3]):
		print(item)
