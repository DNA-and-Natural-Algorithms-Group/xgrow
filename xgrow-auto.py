#!/bin/env python

from numpy import *
import os

def run_nowait(ts,pl):
	temp = os.tempnam()
	pid = os.spawnvp(os.P_NOWAIT, 'xgrow', ['xgrow', ts] + pl + ["tracefile=%s" % temp])
	return (pid, temp)

def read_trace(fn):
	return from_file(fn,sep=" ").
