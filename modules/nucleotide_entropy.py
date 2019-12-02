from __future__ import division
from collections import deque
from decimal import Decimal
from math import log
import sys
import itertools

from modules.codon_probability import CodonProbability


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
		self.entropy_rev = [deque([]),deque([]),deque([]),deque([])]

		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		self.amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.translate = dict(zip(codons, self.amino_acids))

		self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-'}

		
		self.codon_probs = CodonProbability(nucleotides)

		for i, nucl in enumerate(nucleotides):
			# find nucleotide frequency
			self.frame = next(self.frames)
			self.add_base(nucl)

			# find shannon entropy
			#se = self.trinucleotide_entropy(self.frequency[self.frame])
			#se = self.dinucleotide_entropy(self.frequency[self.frame])
			#print(i, nucl, self.codon)	
			se = self.peptide_entropy(self.frequency[self.frame])
			self.entropy[self.frame].append(se)

			#se = self.dinucleotide_entropy(self.reverse_frequencies(self.frequency[self.frame]))
			se = self.peptide_entropy(self.reverse_frequencies(self.frequency[self.frame]))
			self.entropy_rev[self.frame].append(se)
			

	def get(self):
		for _ in range(self.window//2):
			for frame in [1,2,3]:
				self.entropy[frame].popleft()
		for _ in range(self.window//2):
			for frame in [1,2,3]:
				self.entropy_rev[frame].popleft()
		return self.entropy, self.entropy_rev

	def reverse_frequencies(self, dictionary):
		new_dict = dict()
		for key in dictionary:
			new_key = tuple([self.complement[t] for t in key[::-1]])
			new_dict[new_key] = dictionary[key]
		return new_dict

	def add_base(self, base):

		# set newest codon
		self.codon.popleft()
		self.codon.append(base)

 		# add newest codon
		self.codons[self.frame].append(tuple(self.codon))
		self.frequency[self.frame][tuple(self.codon)] += 1
		# misc
		#aa = self.translate.get("".join(self.codon), "-")
		#self.aa[self.frame][aa] += 1

 		# pop oldest codon
		popped_codon = self.codons[self.frame].popleft()
		self.frequency[self.frame][popped_codon] -= 1
		# misc
		#aa = self.translate.get("".join(popped_codon), "-")
		#self.aa[self.frame][aa] -= 1

	def peptide_entropy(self, dictionary):
		se = 0;
		#se0 = 0;
		new_dict = dict()
		total = 0
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				aa = self.translate[''.join(key)]
				#se0 += -self.codon_probs.probability(aa) * log(self.codon_probs.probability(aa))
				count = new_dict.get(aa, 0)
				new_dict[aa] = dictionary[key] + count 
				#total += dictionary[key] + count
		for key in new_dict:
			new_dict[key] = new_dict[key] / self.amino_acids.count(key)
			total += new_dict[key]
		for key in new_dict:
			p = new_dict[key] / 22 #total
			se += -p * log(p)
		return se

	def dinucleotide_entropy(self, dictionary):
		se = 0;
		new_dict = dict()
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				count = new_dict.get(key[0:2] , 0)
				new_dict[key[0:2]] = count + 1
		for key in new_dict:
			p = new_dict[key] / 16
			se += -p * log(p)
		return se

	def trinucleotide_entropy(self, dictionary, direction = 'forward'):
		se = 0;
		for key in dictionary:
			if '-' not in key and dictionary[key]:
				p = dictionary[key] / 64
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

