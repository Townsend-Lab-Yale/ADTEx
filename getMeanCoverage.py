#! /usr/bin/python2.7

# ----------------------------------------------------------------------#
# Copyright (c) 2013, Kaushalya Amarasinghe.
#
# > Source License <
# This file is part of ADTEx.
#
#    ADTEx is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ADTEx is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ADTEx.  If not, see <http://www.gnu.org/licenses/>.
#
# 
#-----------------------------------------------------------------------#

import subprocess
import shlex
import os

def getMeanCoverage(infile,outfile):
	inF = open(infile)
	outF = open(outfile,"w")

	tot = 0
	count = 0

	for line in inF:
		l = line.rstrip().split("\t")
		ll = len(l)
		
		tot = tot + float(l[ll-1])
		if (float(l[ll-1]) <10):
			count = count + 1
	
		if (int(l[2]) - int(l[1]) == int(l[ll-2])):
			#if count > 10:
			if count > -1:
				chr = str(l[0])
				start  = str(l[1])
				end = str(l[2])
				mean = round(float(tot)/int(l[ll-2]),0)
				outF.write(chr+"\t"+start+"\t"+end+"\t"+str(mean)+"\n")
			tot = 0
			count =0
				
	inF.close()
	outF.close()

	#END
