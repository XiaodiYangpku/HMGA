#!/usr/bin/env python

target_group={}
r=open('targets_alltype.txt')
for line in r.readlines():
    line=line.strip('\n')
    target=line.split('\t')[0]
    group=line.split('\t')[1]
    if target not in target_group:target_group[target]=[group]
    else:
        if group not in target_group[target]:target_group[target].append(group)

w=open('2_target_groups.txt', 'w')
for target in target_group:
    w.write(target+'\t'+','.join(target_group[target])+'\n')
