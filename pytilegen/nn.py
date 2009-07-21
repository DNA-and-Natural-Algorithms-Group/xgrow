from numpy import *

"""The start of a nearest neighbor sticky end interaction algorithm
implementation similar to that in DNAdesign, but written in a more convoluted
but likely faster manner. This will eventually be used with tilegen.py to allow
for calculation of interactions from sequences."""

nn_dG_matrix = array([-1.02,-1.43,-1.16,-0.73,-1.38,-1.77,-2.09,-1.16,-1.46,-2.28,-1.77,-1.43,-0.60,-1.46,-1.38,-1.02])

brep = {'a': 0, 'c': 1, 'g': 2, 't': 3}
srep = {0: 'a', 1: 'c', 2: 'g', 3: 't'}

class pairseq:
	def __init__(s,seqstring,in_brep=False):
		if not in_brep:
			s.d=array([4*brep[seqstring[x]]+brep[seqstring[x+1]] for x in range(0,len(seqstring)-1)]);
		else:
			s.d=seqstring;
	def __str__(s):
		st = ""
		st += srep[s.d[0]/4]
		for n in mod(s.d,4):
			st += srep[n]
		return st
	def __array__(s):
		return s.d
	def __repr__(s):
		return "pairseq('%s')" % s.__str__()
	def comp(s):
		return pairseq(invert(s.d)%16,in_brep=True)
	def revcomp(s):
		s2 = flipud(invert(s.d)%16)
		s2 = s2/4 + 4*(s2%4)
		return pairseq(s2,in_brep=True)
	def reverse(s):
		s2 = flipud(s.d)
		return pairseq(s2/4 + 4*(s2%4),in_brep=True)
	def __add__(a,b):
		if a.__class__ == b.__class__ == pairseq:
			return pairseq(concatenate((a.d,array([4*(a.d[-1]%4) + b.d[0]/4]),b.d)),in_brep=True)
		elif a.__class__ == pairseq and b.__class__ == str:
			return pairseq(concatenate((a.d,array([4*(a.d[-1]%4) + brep[b]]))),in_brep=True)
	def __radd__(b,a):
		if a.__class__ == str and b.__class__ == pairseq:
			return pairseq(concatenate((array([4*brep[a] + b.d[0]/4]),b.d)),in_brep=True)
	to_str = __str__

