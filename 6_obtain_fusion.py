#!/usr/bin/env python

import os, sys
import numpy as np

def get_sampleid_infor(clinicalfile):
    r=open(clinicalfile)
    sampleid_infor={}
    r.readline()
    for line in r.readlines():
        line=line.strip('\n')
        infor=line.split('\t')
        sampleid=infor[1]
        sampleid_infor[sampleid]=infor
    return(sampleid_infor)

def get_genome_infor(genomefile):
    r=open(genomefile)
    panel_bp={}
    r.readline()
    for line in r.readlines():
        line=line.strip('\n')
        infor=line.split('\t')
        chr=infor[0]
        chr_s=infor[1]
        chr_e=infor[2]
        gene=infor[3]
        panel=infor[5]
        length=int(chr_e)-int(chr_s)+1
        if panel not in panel_bp:panel_bp[panel]=length
        else:panel_bp[panel]+=length
    r.close()
    return(panel_bp)
        
def get_cancer_infor(dividefile):
    r=open(dividefile)
    r.readline()
    cancer_type={}
    for line in r.readlines():
        line=line.strip()
        if len(line.split())==1:continue
        if line.split('\t')[0]=='MS':continue
        try:
            ctype=line.split('\t')[0]
            detail=line.split('\t')[-1]
            cancer_type[detail]=ctype
        except:
            print(line)
    r.close()
    return(cancer_type)


def get_sampleid_muts(fusionfile, sampleid_infor):
    r=open(fusionfile)
    sampleids=r.readline().strip('\n').split('\t')[1:]
    sampleid_muts_all={}
    panel_fusionlist={}
    for line in r.readlines():
        line=line.strip('\n')
        infor=line.split('\t')
        gene=infor[0]
        sampleid=infor[3]
        fusionoriginal=infor[4]
        fusions=fusionoriginal.split()[0]
        fusion=fusions
        if '-' in fusions:
            try:
                if '-'+gene in fusions:
                    a=fusions.split(gene)[0].strip('-').strip()
                else:
                    a=fusions.split(gene)[1].strip('-').strip()
            except:
                #print(gene,fusions,'aaaa')
            fusion='-'.join(sorted([a,gene]))
            #print(fusion,'fu',line)
            if '-intragenic' in fusions:fusion=gene+'-intragenic'
            if '-intergenic' in fusions:fusion=gene+'-intergenic'
            if 'rearrangement' in fusions:fusion=fusion+' rearrangement'
            elif 'truncation' in fusions:fusion=fusion+' truncation'
            elif 'deletion' in fusions:fusion=fusion+' deletion'
            elif 'duplication' in fusions:fusion=fusion+' duplication'
            elif 'inversion' in fusions:fusion=fusion+' inversion'
            elif 'structural variant' in fusions:fusion=fusion+' structural variant'
            elif 'rearranged' in fusions:fusion=fusion
            elif 'fusion' in fusions:fusion=fusion
            #if 'CIITAin'in fusion:print(line,'cuo',a)

        try:panel=sampleid_infor[sampleid][5]
        except:continue
        if panel not in panel_fusionlist:panel_fusionlist[panel]=[fusion]
        else:
            if fusion not in panel_fusionlist[panel]:panel_fusionlist[panel].append(fusion)
        if sampleid not in sampleid_muts_all:
            sampleid_muts_all[sampleid]=[fusion]
        else:
            if fusion not in sampleid_muts_all[sampleid]:sampleid_muts_all[sampleid].append(fusion)

    r.close()
    return(sampleid_muts_all, panel_fusionlist)
def get_cancer_type(cancertypefile):
    r=open(cancertypefile)
    c_total={}
    for line in r.readlines():
        line=line.strip()
        c_total[line.split('\t')[1]]=line.split('\t')[0]
    return(c_total)

def get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, panel_fusionlist, c_total, outputfile):
    cancer_oa_mut, cancer_cya_mut = {}, {}
    cancer_gene_oa_samplenumber={}
    cancer_gene_cya_samplenumber={}
    ctypes=[]
    for sampleid in sampleid_infor:
        panel=sampleid_infor[sampleid][5]
        cancer=sampleid_infor[sampleid][7]
        if panel not in panel_fusionlist:continue
        try:ctype=cancer_type[cancer]
        except:
#            print(cancer)
            continue
        if ctype not in ctypes:ctypes.append(ctype)
        age=sampleid_infor[sampleid][2]
        if '<' in age:age=age.strip('<')
        if '>' in age:age=age.strip('>')
        if ctype not in cancer_gene_oa_samplenumber:cancer_gene_oa_samplenumber[ctype]={}
        if ctype not in cancer_gene_cya_samplenumber:cancer_gene_cya_samplenumber[ctype]={}
        if ctype not in cancer_oa_mut:cancer_oa_mut[ctype]={}
        if ctype not in cancer_cya_mut:cancer_cya_mut[ctype]={}
        for gene in panel_fusionlist[panel]:
            if gene not in cancer_gene_cya_samplenumber[ctype]:cancer_gene_cya_samplenumber[ctype][gene]=0
            if gene not in cancer_gene_oa_samplenumber[ctype]:cancer_gene_oa_samplenumber[ctype][gene]=0
            if int(age)>40:
                cancer_gene_oa_samplenumber[ctype][gene]+=1
            else:
                cancer_gene_cya_samplenumber[ctype][gene]+=1
            if gene not in cancer_oa_mut[ctype]:cancer_oa_mut[ctype][gene]=0
            if gene not in cancer_cya_mut[ctype]:cancer_cya_mut[ctype][gene]=0
        if int(age)>40:
            if sampleid in sampleid_muts_all:
                for fusion in sampleid_muts_all[sampleid]:
                    try:
                        cancer_oa_mut[ctype][fusion]+=1
                    except:
                        continue
        else:
            if sampleid in sampleid_muts_all:
                for fusion in sampleid_muts_all[sampleid]:
                    try:
                        cancer_cya_mut[ctype][fusion]+=1
                    except:
                        continue
    w=open(outputfile, 'w')
    w.write('cancer\tsubtype\tgene\tCYA_mut\tCYA_wt\tOA_mut\tOA_wt\tCYA_mut_ratio\tOA_mut_ratio\tDiff_CYA_OA\tHigher_in_CYA\tP\n')
    for ctype in ctypes:
        for gene in cancer_cya_mut[ctype]:
            cya_mut=cancer_cya_mut[ctype][gene]
            cya_total=cancer_gene_cya_samplenumber[ctype][gene]
            cya_wt=cya_total-cya_mut
            oa_mut=cancer_oa_mut[ctype][gene]
            oa_total=cancer_gene_oa_samplenumber[ctype][gene]
            oa_wt=oa_total-oa_mut
            try:
                cya_ratio='%.6f'%(cya_mut/cya_total)
            except:
                cya_ratio='NA'
            try:
                oa_ratio='%.6f'%(oa_mut/oa_total)
            except:
                oa_ratio='NA'
            try:
                diff=cya_mut/cya_total-oa_mut/oa_total
            except:
                diff='NA'
            if diff=='NA':cya_status='NA'
            else:
                if diff>0:
                    cya_status='1'
                    diff='%.6f'%abs(diff)
                else:
                    cya_status='0'
                    diff='%.6f'%abs(diff)
            w.write(c_total[ctype]+'\t'+ctype+'\t'+gene+'\t'+str(cya_mut)+'\t'+str(cya_wt)+'\t'+str(oa_mut)+'\t'+str(oa_wt)+'\t'+cya_ratio+'\t'+oa_ratio+'\t'+diff+'\t'+cya_status+'\t')
            w.write('fisher.test(matrix(c('+str(cya_mut)+','+str(cya_wt)+','+str(oa_mut)+','+str(oa_wt)+'),nrow=2))$p.value\n')
    w.close()


def get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, panel_fusionlist, c_total, outputfile):
    cancer_oa_mut, cancer_cya_mut = {}, {}
    ctypes=[]
    cancer_gene_oa_samplenumber={}
    cancer_gene_cya_samplenumber={}
    #print(panel_fusionlist)
    for sampleid in sampleid_infor:
        panel=sampleid_infor[sampleid][5]
        cancer=sampleid_infor[sampleid][7]
        try:ctype=cancer_type[cancer]
        except:
#            print(cancer)
            continue
        if panel not in panel_fusionlist:
            #print('genepanel without ',panel)
            continue
        ctype=c_total[ctype]
        if ctype not in ctypes:ctypes.append(ctype)
        age=sampleid_infor[sampleid][2]
        panel=sampleid_infor[sampleid][5]
        if ctype not in cancer_gene_oa_samplenumber:cancer_gene_oa_samplenumber[ctype]={}
        if ctype not in cancer_gene_cya_samplenumber:cancer_gene_cya_samplenumber[ctype]={}
        if ctype not in cancer_oa_mut:cancer_oa_mut[ctype]={}
        if ctype not in cancer_cya_mut:cancer_cya_mut[ctype]={}
        if '<' in age:age=age.strip('<')
        if '>' in age:age=age.strip('>')
        for gene in panel_fusionlist[panel]:
            #print(gene)
            if gene not in cancer_gene_cya_samplenumber[ctype]:cancer_gene_cya_samplenumber[ctype][gene]=0
            if gene not in cancer_gene_oa_samplenumber[ctype]:cancer_gene_oa_samplenumber[ctype][gene]=0
            if int(age)>40:
                cancer_gene_oa_samplenumber[ctype][gene]+=1
            else:
                cancer_gene_cya_samplenumber[ctype][gene]+=1
            if gene not in cancer_oa_mut[ctype]:cancer_oa_mut[ctype][gene]=0
            if gene not in cancer_cya_mut[ctype]:cancer_cya_mut[ctype][gene]=0
        if int(age)>40:
            if sampleid in sampleid_muts_all:
                for fusion in sampleid_muts_all[sampleid]:
                    try:
                        cancer_oa_mut[ctype][fusion]+=1
                    except:
                        #print(fusion)
                        continue
        else:
            if sampleid in sampleid_muts_all:
                for fusion in sampleid_muts_all[sampleid]:
                    try:
                        cancer_cya_mut[ctype][fusion]+=1
                    except:
                        #print(fusion)
                        continue
    w=open(outputfile, 'w')
    w.write('cancer\tgene\tCYA_mut\tCYA_wt\tOA_mut\tOA_wt\tCYA_mut_ratio\tOA_mut_ratio\tDiff_CYA_OA\tHigher_in_CYA\tP\n')
    for ctype in ctypes:
        nan=''
        for gene in cancer_cya_mut[ctype]:
            cya_mut=cancer_cya_mut[ctype][gene]
            cya_total=cancer_gene_cya_samplenumber[ctype][gene]
            cya_wt=cya_total-cya_mut
            oa_mut=cancer_oa_mut[ctype][gene]
            oa_total=cancer_gene_oa_samplenumber[ctype][gene]
            oa_wt=oa_total-oa_mut
            try:
                cya_ratio='%.6f'%(cya_mut/cya_total)
            except:
                cya_ratio='NA'
            try:
                oa_ratio='%.6f'%(oa_mut/oa_total)
            except:
                oa_ratio='NA'
            try:
                diff=cya_mut/cya_total-oa_mut/oa_total
            except:
                diff='NA'
            if diff=='NA':
                cya_status='NA'
            else:
                if diff>0:
                    cya_status='1'
                    diff='%.6f'%abs(diff)
                else:
                    cya_status='0'
                    diff='%.6f'%abs(diff)
            w.write(ctype+'\t'+gene+'\t'+str(cya_mut)+'\t'+str(cya_wt)+'\t'+str(oa_mut)+'\t'+str(oa_wt)+'\t'+cya_ratio+'\t'+oa_ratio+'\t'+diff+'\t'+cya_status+'\t')
            w.write('fisher.test(matrix(c('+str(cya_mut)+','+str(cya_wt)+','+str(oa_mut)+','+str(oa_wt)+'),nrow=2))$p.value\n')
    w.close()


def main():
    clinicalfile='./patient_information.txt'
    fusionfile='genie-data/data_fusions_revised.txt'
    genomefile='genie-data/genomic_information.txt'
    dividefile='hematologic_divided.txt'
    outputfile1='results/6_fusion_csubtype.txt'
    outputfile2='results/6_fusion_ctype.txt'
    cancertypefile='results/cancer_type.txt'
    panelfiledir='genie-data/gene_panels/'

    cancer_type=get_cancer_infor(dividefile)
    sampleid_infor=get_sampleid_infor(clinicalfile)
    panel_bp=get_genome_infor(genomefile)
    sampleid_muts_all, panel_fusionlist=get_sampleid_muts(fusionfile, sampleid_infor)
    c_total=get_cancer_type(cancertypefile)
    get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, panel_fusionlist, c_total, outputfile2)
    get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, panel_fusionlist, c_total, outputfile1)

if __name__=='__main__':
    main()
