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


def get_sampleid_muts(cnafile):
    r=open(cnafile)
    sampleids=r.readline().strip('\n').split('\t')[1:]
    sampleid_muts_pos={}
    sampleid_muts_neg={}
    sampleid_muts_all={}
    for line in r.readlines():
        line=line.strip('\n')
        infor=line.split('\t')
        gene=infor[0]
        values=infor[1:]
        for i in range(len(values)):
            value=values[i]
            sampleid=sampleids[i]
            if value=='NA':continue
            if sampleid not in sampleid_muts_all:
                sampleid_muts_all[sampleid]=[gene]
            else:
                if gene not in sampleid_muts_all[sampleid]:sampleid_muts_all[sampleid].append(gene)
            if float(value)>0:
                if sampleid not in sampleid_muts_pos:
                    sampleid_muts_pos[sampleid]=[gene]
                else:
                    if gene not in sampleid_muts_pos[sampleid]:sampleid_muts_pos[sampleid].append(gene)
            elif float(value)<0:
                if sampleid not in sampleid_muts_neg:
                    sampleid_muts_neg[sampleid]=[gene]
                else:
                    if gene not in sampleid_muts_neg[sampleid]:sampleid_muts_neg[sampleid].append(gene)

    r.close()
    return(sampleid_muts_all, sampleid_muts_pos, sampleid_muts_neg)
def get_cancer_type(cancertypefile):
    r=open(cancertypefile)
    c_total={}
    for line in r.readlines():
        line=line.strip()
        c_total[line.split('\t')[1]]=line.split('\t')[0]
    return(c_total)

def get_genepanel_genes(panelfiledir):
    panellist=os.listdir(panelfiledir)
    panel_genelist={}
    for panel in panellist:
        r=open(panelfiledir+panel)
        genelist=r.readlines()[2].strip('\n').split('\t')[1:]
        r.close()
        panel_genelist[panel.split('_')[-1].strip('.txt')]=genelist
    return(panel_genelist)

def get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, sampleid_muts_pos, sampleid_muts_neg, c_total, outputfile, panel_genelist):
    cancer_oa_mut, cancer_cya_mut = {}, {}
    cancer_gene_oa_samplenumber={}
    cancer_gene_cya_samplenumber={}
    ctypes=[]
    for sampleid in sampleid_infor:
        panel=sampleid_infor[sampleid][5]
        cancer=sampleid_infor[sampleid][7]
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
        for gene in panel_genelist[panel]:
            if gene not in cancer_gene_cya_samplenumber[ctype]:cancer_gene_cya_samplenumber[ctype][gene]=0
            if gene not in cancer_gene_oa_samplenumber[ctype]:cancer_gene_oa_samplenumber[ctype][gene]=0
            if int(age)>40:
                cancer_gene_oa_samplenumber[ctype][gene]+=1
            else:
                cancer_gene_cya_samplenumber[ctype][gene]+=1
            if gene not in cancer_oa_mut[ctype]:cancer_oa_mut[ctype][gene]=0
            if gene not in cancer_cya_mut[ctype]:cancer_cya_mut[ctype][gene]=0
        if int(age)>40:
            if sampleid in sampleid_muts_neg:
                for gene in sampleid_muts_neg[sampleid]:
                    try:
                        cancer_oa_mut[ctype][gene]+=1
                    except:
                        continue
        else:
            if sampleid in sampleid_muts_neg:
                for gene in sampleid_muts_neg[sampleid]:
                    try:
                        cancer_cya_mut[ctype][gene]+=1
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
#            w.write(c_total[ctype]+'\t'+ctype+'\t'+gene+'\t'+str(cya_mut)+'\t'+str(cya_wt)+'\t'+str(oa_mut)+'\t'+str(oa_wt)+'\t')
            w.write('fisher.test(matrix(c('+str(cya_mut)+','+str(cya_wt)+','+str(oa_mut)+','+str(oa_wt)+'),nrow=2))$p.value\n')
    w.close()


def get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, sampleid_muts_pos, sampleid_muts_neg, c_total, outputfile, panel_genelist):
    cancer_oa_mut, cancer_cya_mut = {}, {}
    ctypes=[]
    cancer_gene_oa_samplenumber={}
    cancer_gene_cya_samplenumber={}
    #print(panel_genelist)
    for sampleid in sampleid_infor:
        panel=sampleid_infor[sampleid][5]
        cancer=sampleid_infor[sampleid][7]
        try:ctype=cancer_type[cancer]
        except:
#            print(cancer)
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
        if panel not in panel_genelist:
            #print('genepanel without ',panel)
            continue
        for gene in panel_genelist[panel]:
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
            if sampleid in sampleid_muts_neg:
                for gene in sampleid_muts_neg[sampleid]:
                    try:
                        cancer_oa_mut[ctype][gene]+=1
                    except:
                        #print(gene)
                        continue
        else:
            if sampleid in sampleid_muts_neg:
                for gene in sampleid_muts_neg[sampleid]:
                    try:
                        cancer_cya_mut[ctype][gene]+=1
                    except:
                        #print(gene)
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
    clinicalfile='./patient_information_have_dod_or_contact.txt'
    cnafile='genie-data/data_CNA.txt'
    genomefile='genie-data/genomic_information.txt'
    dividefile='hematologic_divided.txt'
    outputfile1='results/5_cna_neg_csubtype.txt'
    outputfile2='results/5_cna_neg_ctype.txt'
    cancertypefile='results/cancer_type.txt'
    panelfiledir='genie-data/gene_panels/'

    panel_genelist=get_genepanel_genes(panelfiledir)
    cancer_type=get_cancer_infor(dividefile)
    sampleid_infor=get_sampleid_infor(clinicalfile)
    panel_bp=get_genome_infor(genomefile)
    sampleid_muts_all, sampleid_muts_pos, sampleid_muts_neg=get_sampleid_muts(cnafile)
    c_total=get_cancer_type(cancertypefile)
    get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, sampleid_muts_pos, sampleid_muts_neg, c_total, outputfile2, panel_genelist)
    get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts_all, sampleid_muts_pos, sampleid_muts_neg, c_total, outputfile1, panel_genelist)

if __name__=='__main__':
    main()
