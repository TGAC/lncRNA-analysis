#!/usr/bin/env python

import sys
attr = {}
for line in sys.stdin:
	if(line.split('\t')[2] == "transcript"):
		for item in line.split('\t')[8].strip('\n').split('; '):
			(key,value) = item.split(' ')
			attr[key] = value
#			print(attr)
		sys.stdout.write("%s\n" % (attr['gene_id'].strip('\"') + '\t' + attr['transcript_id'].strip('\"')))
		#print(attr['gene_id'])
		#print(attr['transcript_id'])

#    attr = dict(item.strip().split(' ') for item in line.split('\t')[8].strip('\n').split(' ') if item)
 #   sys.stdout.write("%s\n" % (attr['gene_id'].strip('\"') + '\t' + attr['transcript_id'].strip('\"')))

