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


def get_sampleid_muts(mutfile):
    r=open(mutfile)
    r.readline()
    sampleid_muts={}
    for line in r.readlines():
        line=line.strip('\n')
        infor=line.split('\t')
        gene=infor[0]
        chr=infor[4]
        chr_s=infor[5]
        chr_e=infor[6]
        mut_type1=infor[8]
        mut_type2=infor[9]
        var_type=infor[10]
        ref=infor[11]
        alt=infor[13]
        dbsnp=infor[14]
        sampleid=infor[16]
        t_ref_count=infor[33]
        t_alt_count=infor[34]
        codon=infor[37] #ENST00000256078.4:c.34G>T 
        protein=infor[38] # p.Gly12Cys
        pro_short=infor[39] # `p.G12C
        transcriptome_id=infor[40] #ENST00000256078
        refseq=infor[41] # NM_
        pro_position=infor[42] # 12
        codons=infor[43] # Ggt/Tgt
        exon_number=infor[44] # 15/18   
        gnomAD_AF=infor[45] #3.97994e-06
        polyphen=infor[55] #Polyphen_Prediction     Polyphen_Score  SIFT_Prediction SIFT_Score
        sift=infor[57] #probably_damaging       0.991   deleterious     0.04 
        t_depth=infor[61]
        if sampleid not in sampleid_muts:
            sampleid_muts[sampleid]=[gene]
        else:
            if gene not in sampleid_muts[sampleid]:sampleid_muts[sampleid].append(gene)
    r.close()
    return(sampleid_muts)
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

def get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile, panel_genelist):
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
            if sampleid in sampleid_muts:
                for gene in sampleid_muts[sampleid]:
                    try:
                        cancer_oa_mut[ctype][gene]+=1
                    except:
                        continue
        else:
            if sampleid in sampleid_muts:
                for gene in sampleid_muts[sampleid]:
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


def get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile, panel_genelist):
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
            if sampleid in sampleid_muts:
                for gene in sampleid_muts[sampleid]:
                    try:
                        cancer_oa_mut[ctype][gene]+=1
                    except:
                        #print(gene)
                        continue
        else:
            if sampleid in sampleid_muts:
                for gene in sampleid_muts[sampleid]:
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
    clinicalfile='./patient_information.txt'
    mutfile='genie-data/data_mutations_extended.txt'
    genomefile='genie-data/genomic_information.txt'
    dividefile='hematologic_divided.txt'
    outputfile1='results/4_mutations_csubtype.txt'
    outputfile2='results/4_mutations_ctype.txt'
    cancertypefile='results/cancer_type.txt'
    panelfiledir='genie-data/gene_panels/'

    panel_genelist=get_genepanel_genes(panelfiledir)
    cancer_type=get_cancer_infor(dividefile)
    sampleid_infor=get_sampleid_infor(clinicalfile)
    panel_bp=get_genome_infor(genomefile)
    sampleid_muts=get_sampleid_muts(mutfile)
    c_total=get_cancer_type(cancertypefile)
    get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile2, panel_genelist)
    get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile1, panel_genelist)

if __name__=='__main__':
    main()
