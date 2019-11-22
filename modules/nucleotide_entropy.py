from __future__ import division
from collections import deque
from decimal import Decimal
from math import log
import sys
import itertools

class Base:
    def __init__(self, nucl):
       self.nucl = nucl
       self.next = None
 
class Codon:
    def __init__(self):
        self.head = None
        self.tail = None
 
    def push(self,  nucl):
        if self.head is None:
            self.head = Base(nucl)
        else:
            new_nucl = Base(nucl)
            new_nucl.next = self.head
            self.head = new_nucl

    def pull(self):
        pass
 
    def pop(self):
        if self.head is None:
            return None
        else:
            popped = self.head.nucl
            self.head = self.head.next
            return popped



class NucleotideEntropy:
	def __init__(self, nucleotides, window = 150):
		self.window = window//3
		self.frames = itertools.cycle([1, 2, 3])
		self.frame = None

		self.codon = deque(['-','-','-'])

		self.codons = [None] * 4
		self.codons[1] = deque([tuple(['-','-','-'])] * self.window)
		self.codons[2] = deque([tuple(['-','-','-'])] * self.window)
		self.codons[3] = deque([tuple(['-','-','-'])] * self.window)		

		self.frequency = [None] * 4
		self.frequency[1] = self._init_dict()
		self.frequency[2] = self._init_dict()
		self.frequency[3] = self._init_dict()

		self.entropy = [deque([]),deque([]),deque([]),deque([])]

		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.translate = dict(zip(codons, amino_acids))

		self.aa = [None] * 4
		self.aa[1] = self._init_aadict()
		self.aa[2] = self._init_aadict()
		self.aa[3] = self._init_aadict()

		for nucl in nucleotides:
			# find nucleotide frequency
			self.frame = next(self.frames)
			self.add_base(nucl)

			# find shannon entropy
			#se = self.shannon_entropy(self.frequency[self.frame])
			se = self.shannon_entropy2(self.aa[self.frame])
			self.entropy[self.frame].append(se)

	def entropy_at(self):
		return self.entropy

	def add_base(self, base):

		# set newest codon
		self.codon.popleft()
		self.codon.append(base)

 		# add newest codon
		self.codons[self.frame].append(tuple(self.codon))
		self.frequency[self.frame][tuple(self.codon)] += 1
		# misc
		aa = self.translate.get("".join(self.codon), "-")
		self.aa[self.frame][aa] += 1

 		# pop oldest codon
		popped_codon = self.codons[self.frame].popleft()
		self.frequency[self.frame][popped_codon] -= 1
		# misc
		aa = self.translate.get("".join(popped_codon), "-")
		self.aa[self.frame][aa] -= 1

	def translate(self):
		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		codon_table = dict(zip(codons, amino_acids))

	def shannon_entropy(self, dictionary):
		se = 0;
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				p = dictionary[key] / 64
				se += -p * log(p)
		return se


	def shannon_entropy2(self, dictionary):
		se = 0;
		for key in dictionary:
			if key != '-' and dictionary[key]:
				p = dictionary[key] / self.window
				se += -p * log(p)
		return se

	def _init_dict(self):
		freq_dict = dict()
		for first in 'ATGC-':
			for second in 'ATGC-':
				for third in 'ATGC-':
					freq_dict[(tuple([first, second, third]))] = 0
		return freq_dict
				
	def _init_aadict(self):
		freq_dict = dict()
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-'
		for aa in amino_acids:
			freq_dict[aa] = 0
		return freq_dict

	def _close(self):
		for _ in range(self.window//2):
			for frame in [1,2,3]:
				self.total[frame].popleft()
				item = self.bases[frame].popleft()
				self.frequency[frame][item] -= 1
				self.total[frame].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])

	def get(self):
		self._close()
		self.freq.append( [20,20,20] )
		for i in range(len(self.total[3])-1):
			self.freq.append( [self.total[1][i],self.total[2][i],self.total[3][i]] )
			self.freq.append( [self.total[2][i],self.total[3][i],self.total[1][i+1]] )
			self.freq.append( [self.total[3][i],self.total[1][i+1],self.total[2][i+1]] )
		i += 1
		self.freq.append( [self.total[1][i],self.total[2][i],self.total[3][i]] )
		if(i < len(self.total[1])-1):
			self.freq.append( [self.total[2][i],self.total[3][i],self.total[1][i+1]] )
		if(i < len(self.total[2])-1):
			self.freq.append( [self.total[3][i],self.total[1][i+1],self.total[2][i+1]] )
		return self.freq

