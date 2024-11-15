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
            sampleid_muts[sampleid].append(gene)
    r.close()
    return(sampleid_muts)
def get_cancer_type(cancertypefile):
    r=open(cancertypefile)
    c_total={}
    for line in r.readlines():
        line=line.strip()
        c_total[line.split('\t')[1]]=line.split('\t')[0]
    return(c_total)

def get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile):
    cancer_oa_tmb, cancer_cya_tmb = {}, {}
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
        if sampleid in sampleid_muts:
            mutnum=len(sampleid_muts[sampleid])
        else:
            mutnum=0
        tmb=mutnum/(panel_bp[panel]/1000000)
        if panel_bp[panel]/1000000 <0.9:continue
        #print(cancer+'\t'+ctype,sampleid,mutnum,mutnum/(panel_bp[panel]/1000000))
        if '<' in age:age=age.strip('<')
        if '>' in age:age=age.strip('>')
        if int(age)>39:
            if ctype not in cancer_oa_tmb:cancer_oa_tmb[ctype]=[tmb]
            else:cancer_oa_tmb[ctype].append(tmb)
        else:
            if ctype not in cancer_cya_tmb:cancer_cya_tmb[ctype]=[tmb]
            else:cancer_cya_tmb[ctype].append(tmb)
    w=open(outputfile, 'w')
    for ctype in ctypes:
        nan=''
        w.write(c_total[ctype]+'\t'+ctype+'\t')
        if ctype in cancer_cya_tmb:
            average='%.6f'%(sum(cancer_cya_tmb[ctype])/len(cancer_cya_tmb[ctype]))
            median=np.median(cancer_cya_tmb[ctype])
            cancer_cya_tmb[ctype]=['%.6f'%(i) for i in cancer_cya_tmb[ctype]]
            w.write(str(median)+'\t'+average+'\t')
        else:
            w.write('NA\tNA\t')
            nan='NA'
        if ctype in cancer_oa_tmb:
            average='%.6f'%(sum(cancer_oa_tmb[ctype])/len(cancer_oa_tmb[ctype]))
            median=np.median(cancer_oa_tmb[ctype])
            cancer_oa_tmb[ctype]=['%.6f'%(i) for i in cancer_oa_tmb[ctype]]
            w.write(str(median)+'\t'+average+'\t')
        else:
            nan='NA'
            w.write('NA\tNA\t')
        if nan=='NA':w.write('NA\n')
        else:w.write('wilcox.test(c('+','.join(cancer_cya_tmb[ctype])+'),c('+','.join(cancer_oa_tmb[ctype])+'))$p.value\n')
    w.close()


def get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile):
    cancer_oa_tmb, cancer_cya_tmb = {}, {}
    ctypes=[]
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
        if sampleid in sampleid_muts:
            mutnum=len(sampleid_muts[sampleid])
        else:
            mutnum=0
        tmb=mutnum/(panel_bp[panel]/1000000)
        if panel_bp[panel]/1000000 <0.9:continue
        #print(cancer+'\t'+ctype,sampleid,mutnum,mutnum/(panel_bp[panel]/1000000))
        if '<' in age:age=age.strip('<')
        if '>' in age:age=age.strip('>')
        if int(age)>39:
            if ctype not in cancer_oa_tmb:cancer_oa_tmb[ctype]=[tmb]
            else:cancer_oa_tmb[ctype].append(tmb)
        else:
            if ctype not in cancer_cya_tmb:cancer_cya_tmb[ctype]=[tmb]
            else:cancer_cya_tmb[ctype].append(tmb)
    w=open(outputfile, 'w')
    for ctype in ctypes:
        nan=''
        w.write(ctype+'\t')
        if ctype in cancer_cya_tmb:
            average='%.6f'%(sum(cancer_cya_tmb[ctype])/len(cancer_cya_tmb[ctype]))
            median=np.median(cancer_cya_tmb[ctype])
            cancer_cya_tmb[ctype]=['%.6f'%(i) for i in cancer_cya_tmb[ctype]]
            w.write(str(median)+'\t'+average+'\t')
        else:
            nan='NA'
            w.write('NA\tNA\t')
        if ctype in cancer_oa_tmb:
            average='%.6f'%(sum(cancer_oa_tmb[ctype])/len(cancer_oa_tmb[ctype]))
            median=np.median(cancer_oa_tmb[ctype])
            cancer_oa_tmb[ctype]=['%.6f'%(i) for i in cancer_oa_tmb[ctype]]
            w.write(str(median)+'\t'+average+'\t')
        else:
            nan='NA'
            w.write('NA\tNA\n')
        if nan=='NA':w.write('NA\t')
        else:w.write('wilcox.test(c('+','.join(cancer_cya_tmb[ctype])+'),c('+','.join(cancer_oa_tmb[ctype])+'))$p.value\n')
    w.close()


def main():
    clinicalfile='./patient_information.txt'
    mutfile='genie-data/data_mutations_extended.txt'
    genomefile='genomic_information.txt'
    dividefile='hematologic_divided.txt'
    outputfile1='results/3_TMB_csubtype.txt'
    outputfile2='results/3_TMB_ctype.txt'
    cancertypefile='results/cancer_type.txt'

    cancer_type=get_cancer_infor(dividefile)
    sampleid_infor=get_sampleid_infor(clinicalfile)
    panel_bp=get_genome_infor(genomefile)
    sampleid_muts=get_sampleid_muts(mutfile)
    c_total=get_cancer_type(cancertypefile)
    get_oa_cya_tmb_ctypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile2)
    get_oa_cya_tmb_csubtypes(cancer_type, sampleid_infor, panel_bp, sampleid_muts, c_total, outputfile1)

if __name__=='__main__':
    main()
