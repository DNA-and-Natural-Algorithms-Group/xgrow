#!/bin/env python

from numpy import *
import os

def run_nowait(ts,pl):
	temp = os.tempnam()
	pid = os.spawnvp(os.P_NOWAIT, 'xgrow', ['xgrow', ts] + pl + ["tracefile=%s" % temp])
	return (pid, temp)

def read_trace(fn):
	return from_file(fn,sep=" ").reshape((-1,10))
	# Remember that trace is [Gmc,Gse,ratek,t,tiles,mismatches,events,perimeter,G,dG_bonds,text?]

def run_xgrows(tileset,pll,conc):
	procdict = {}
	reslist = []
	pll.reverse()
	nprocs = 0
	i = 0
	while nprocs < conc:
		(pid, temp) = run_nowait(tileset,pl)
		procdict += {pid: (i,temp)}
		i += 1
		nprocs += 1
	while pll:
		fpid = os.wait()[0]
		pi = procdict[fpid]
		reslist += [(pi[0],read_trace(pi[1]))]
		(pid, temp) = run_nowait(tileset,pl)
		procdict += {pid: (i,temp)}
		i+=1
	while nprocs > 0:
		fpid = os.wait()[0]
		pi = procdict[fpid]
		reslist += [(pi[0],read_trace(pi[1]))]
		nprocs -= 1
	reslist.sort()
	return reslist

def gen_params_gmc_gse(GseVals,GmcVals,otherparams):
	return [ ['Gmc=%d' % Gse,'Gse=%d' % Gmc]+otherparams for Gse in GseVals for Gmc in GmcVals ]
