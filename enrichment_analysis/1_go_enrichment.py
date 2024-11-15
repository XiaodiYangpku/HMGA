#!/usr/bin/env python
#-*-coding:utf-8 -*-

# GO enrichment

#step1 get module gene list
def get_module_targets(modulefile):
    r=open(modulefile)
    r.readline()
    genelist=[]
    for line in r.readlines():
        line=line.strip()
        human_gene=line
        if human_gene not in genelist:
            genelist.append(human_gene)
    return(genelist)

#step2 get human->GOCC,GOBP,GOMF
def get_human_GO(gofile):
    r=open(gofile)
    human_gocc,human_gobp,human_gomf={},{},{}
    Ngo_num={}
    for line in r.readlines():
        line=line.strip()
        if line[0]=='!':continue
        humanuniprotid=line.split('\t')[0]
        goid=line.split('\t')[1]
        gotype=line.split('\t')[2]
        if goid not in Ngo_num:
            Ngo_num[goid]=1
        else:
            Ngo_num[goid]+=1
        if gotype=='C':
            if humanuniprotid not in human_gocc:
                human_gocc[humanuniprotid]=[goid]
            else:
                if goid not in human_gocc[humanuniprotid]:
                    human_gocc[humanuniprotid].append(goid)
        elif gotype=='F':
            if humanuniprotid not in human_gomf:
                human_gomf[humanuniprotid]=[goid]
            else:
                if goid not in human_gomf[humanuniprotid]:
                    human_gomf[humanuniprotid].append(goid)
        else:
            if humanuniprotid not in human_gobp:
                human_gobp[humanuniprotid]=[goid]
            else:
                if goid not in human_gobp[humanuniprotid]:
                    human_gobp[humanuniprotid].append(goid)
    return(human_gobp,human_gomf,human_gocc,Ngo_num)

if __name__=='__main__':
    human_gobp,human_gomf,human_gocc,Ngo_num=get_human_GO('goa_uniq')
    Nbp=len(human_gobp)
    Nmf=len(human_gomf)
    Ncc=len(human_gocc)
    w=open('./1_goenrichment_input','w')
    targetfiles=['../targets/targets_myeloid_CYA0.5.txt','../targets/targets_myeloid_OA0.5.txt','../targets/targets_lymphoid_CYA0.5.txt','../targets/targets_lymphoid_OA0.5.txt']
    for i in targetfiles:
        genelist=get_module_targets(i)
        nbp,ncc,nmf=0,0,0
        humanlist=genelist
        ngobp_num,ngocc_num,ngomf_num={},{},{}
        for eachhuman in humanlist:
            if eachhuman in human_gobp:
                nbp+=1
                for eachgo in human_gobp[eachhuman]:
                    if  eachgo not in ngobp_num:
                        ngobp_num[eachgo]=[eachhuman]
                    else:
                        ngobp_num[eachgo].append(eachhuman)
            if eachhuman in human_gocc:
                ncc+=1
                for eachgo in human_gocc[eachhuman]:
                    if  eachgo not in ngocc_num:
                        ngocc_num[eachgo]=[eachhuman]
                    else:
                        ngocc_num[eachgo].append(eachhuman)
            if eachhuman in human_gomf:
                nmf+=1
                for eachgo in human_gomf[eachhuman]:
                    if  eachgo not in ngomf_num:
                        ngomf_num[eachgo]=[eachhuman]
                    else:
                        ngomf_num[eachgo].append(eachhuman)
        for eachbp in ngobp_num:
            Mbp=Ngo_num[eachbp]
            kbp=len(ngobp_num[eachbp])
            w.write('List'+str(i)+'\tBP\t'+eachbp+'\t'+str(kbp)+'\t'+str(Mbp)+'\t'+str(Nbp)+'\t'+str(nbp)+'\t'+','.join(ngobp_num[eachbp])+'\n')
        for eachcc in ngocc_num:
            Mcc=Ngo_num[eachcc]
            kcc=len(ngocc_num[eachcc])
            w.write('List'+str(i)+'\tCC\t'+eachcc+'\t'+str(kcc)+'\t'+str(Mcc)+'\t'+str(Ncc)+'\t'+str(ncc)+'\t'+','.join(ngocc_num[eachcc])+'\n')
        for eachmf in ngomf_num:
            Mmf=Ngo_num[eachmf]
            kmf=len(ngomf_num[eachmf])
            w.write('List'+str(i)+'\tMF\t'+eachmf+'\t'+str(kmf)+'\t'+str(Mmf)+'\t'+str(Nmf)+'\t'+str(nmf)+'\t'+','.join(ngomf_num[eachmf])+'\n')
