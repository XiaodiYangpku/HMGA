#!/usr/bin/env python

def get_targets(ppifile, genelist, wfile):
    tmp=[]
    r=open(ppifile)
    w=open(wfile, 'w')
    ppis=[]
    targets=[]
    for line in r.readlines():
        line=line.strip('\n')
        pro1=line.split('\t')[0]
        pro2=line.split('\t')[1]
        if pro1==pro2:continue
        if pro1 in genelist:
            if pro2 not in targets:
                targets.append(pro2)
                w.write(pro2+'\n')
                if line not in ppis:ppis.append(line)
        if pro2 in genelist:
            if pro1 not in targets:
                targets.append(pro1)
                w.write(pro1+'\n')
                ppis.append(line)
                if line not in ppis:ppis.append(line)
        tmp.append(line)
        if pro1 in genelist and pro2 in genelist:print(line)
    print('start')
    for t in tmp:
        p1=t.split('\t')[0]
        p2=t.split('\t')[1]
        if p1 in targets or p1 in genelist:
            if p2 in targets or p2 in genelist:
                if t not in ppis:ppis.append(t)
    for ppi in ppis:
        print(ppi)
    r.close()
    w.close()
        
def main():
    ppifile='hppi_genename_score0.5.txt'
    genelist1=['MYD88','KMT2D','BCL2','CD79B','TP53','TET2','CREBBP','DUSP2','DNMT3A','PIM1'] # Lymphoid OA
    genelist2=['NRAS','KRAS','SMARCA4','ID3','PTPN11','RPTOR','MYC'] # Lymphoid CYA
    genelist3=['TET2','SRSF2','ASXL1','DNMT3A','TP53','RUNX1','U2AF1','SF3B1','JAK2','IDH2','ZRSR2'] # Myeloid CA
    genelist4=['KIT','WT1','EP300'] # Myeloid CYA
    genelist=genelist1+genelist2+genelist3+genelist4
    

    wfile='targets_all.txt'
    
    get_targets(ppifile, genelist, wfile)

if __name__=='__main__':
    main()
