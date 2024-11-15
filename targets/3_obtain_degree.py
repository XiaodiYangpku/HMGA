#!/usr/bin/env python

def get_pro_degree(topologyfile):
    r=open(topologyfile)
    r.readline()
    gene_degree={}
    for line in r.readlines():
        line=line.strip()
        gene=line.split()[0]
        degree=line.split()[1]
        betweenness=line.split()[2]
        closeness=line.split()[3]
        gene_degree[gene]=degree
    r.close()
    return(gene_degree)

def get_average_degree(gene_degree, targetfile, outputfile):
    r=open(targetfile)
    group_degree, group_targetnum = {}, {}
    group_eachdegree={}
    for line in r.readlines():
        line=line.strip('\n')
        target=line.split('\t')[0]
        group=line.split('\t')[1]
        if group not in group_degree:
            group_degree[group]=int(gene_degree[target])
            group_eachdegree[group]=[gene_degree[target]]
            group_targetnum[group]=1
            
        else:
            group_degree[group]+=int(gene_degree[target])
            group_targetnum[group]+=1
            group_eachdegree[group].append(gene_degree[target])
    r.close()
    w=open(outputfile, 'w')
    w.write('Group\tDegree\n')
    for group in group_degree:
        for each in group_eachdegree[group]:
            w.write(group+'\t'+each+'\n')
    w.close()
    
def main():
    topologyfile='ppi_topology.txt'
    targetfile='targets_alltype.txt'
    outputfile='3_degree.txt'

    gene_degree=get_pro_degree(topologyfile)
    get_average_degree(gene_degree, targetfile, outputfile)

if __name__=='__main__':
    main()
