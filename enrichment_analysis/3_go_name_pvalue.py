#!/usr/bin/env python
#-*-coding:utf-8 -*-
import os,sys
# get GOid <-> GO name
def get_goid_goname(goobo):
	r=open(goobo)
	gonames=[]
	goid,goname='',''
	goid_goname={}
	filelist=r.readlines()
	for i in range(len(filelist)):
		line=filelist[i]
		if line[0:6]=='id: GO':
			goid=line.strip().split(': ')[1]
			goname=filelist[i+1].strip().split(': ')[1]
			goid_goname[goid]=goname
	r.close()
	return(goid_goname)

if __name__=='__main__':
	goid_goname=get_goid_goname('go.obo')
	r=open('2_GO_pvalue')
	w=open('3_goid_goname_pvalue','w')
	w.write('Virus_entry\tGOtype\tGOid\tGOname\tk\tM\tN\tn\tPvalue\tAdjpvalue\n')
	adjp_other={}
	for line in r.readlines():
		line=line.strip()
		goid=line.split('\t')[2]
		adjp=line.split('\t')[9]
		if len(adjp)>7:
			adjp='%.2f'%float(adjp.split('e')[0])+'e'+adjp.split('e')[1]
		if float(adjp)>0.05:continue
		try:
			if float(adjp) not in adjp_other:
				adjp_other[float(adjp)]=['\t'.join(line.split('\t')[0:3])+'\t'+goid_goname[goid]+'\t'+'\t'.join(line.split('\t')[3:9])+'\t'+adjp]
			else:
				adjp_other[float(adjp)].append('\t'.join(line.split('\t')[0:3])+'\t'+goid_goname[goid]+'\t'+'\t'.join(line.split('\t')[3:9])+'\t'+adjp)
		except:
			print(goid)
	for adjp in sorted(adjp_other):
		for other in adjp_other[adjp]:
			w.write(other+'\n')
