#!/usr/bin/env python

### divided hematological patients into CYA (<=40) and OA (>40) ###
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

def get_cancer_age(inforfile, cancer_type, outputfile, c_total_type):
    r=open(inforfile)
    cancer_cya_oa = {}
    r.readline()
    for line in r.readlines():
        line=line.strip()
        age=line.split('\t')[2]
        cancer=line.split('\t')[7]
        #print(age, cancer)
        if cancer not in cancer_type:
            #print('error', cancer)
            continue
        ctype=cancer_type[cancer]
        if ctype not in cancer_cya_oa:
            cancer_cya_oa[ctype]={}
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']=0
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']=0
            if '>' in age:
                cancer_cya_oa[ctype]['oa']+=1
            elif age == '<18':cancer_cya_oa[ctype]['cya']+=1
            elif int(age)<41:cancer_cya_oa[ctype]['cya']+=1
            else:cancer_cya_oa[ctype]['oa']+=1
        else:
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']=0
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']=0
            if '>' in age:
                cancer_cya_oa[ctype]['oa']+=1
            elif age == '<18':cancer_cya_oa[ctype]['cya']+=1
            elif int(age)<41:cancer_cya_oa[ctype]['cya']+=1
            else:cancer_cya_oa[ctype]['oa']+=1

    r.close()
    w=open(outputfile, 'w')
    w.write('Type\tDisease\tnumber_of_CYA\tnumber_of_OA\tTotal\n')
    for ctype in cancer_cya_oa:
        w.write(c_total_type[ctype]+'\t'+ctype+'\t'+str(cancer_cya_oa[ctype]['cya'])+'\t'+str(cancer_cya_oa[ctype]['oa'])+'\t'+str(cancer_cya_oa[ctype]['cya']+cancer_cya_oa[ctype]['oa'])+'\n')
    w.close()

def get_cancer_sex(inforfile, cancer_type, outputfile, c_total_type):
    r=open(inforfile)
    cancer_cya_oa = {}
    r.readline()
    for line in r.readlines():
        line=line.strip()
        age=line.split('\t')[2]
        cancer=line.split('\t')[7]
        sex=line.split('\t')[9]
        #print(age, cancer)
        if cancer not in cancer_type:
        #    print('error', cancer)
            continue
        ctype=cancer_type[cancer]
        if ctype not in cancer_cya_oa:
            cancer_cya_oa[ctype]={}
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']={}
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']={}
            if '>' in age:
                if sex not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][sex]=1
                else:
                    cancer_cya_oa[ctype]['oa'][sex]+=1
            elif age == '<18':
                if sex not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][sex]=1
                else:
                    cancer_cya_oa[ctype]['cya'][sex]+=1
            elif int(age)<41:
                if sex not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][sex]=1
                else:
                    cancer_cya_oa[ctype]['cya'][sex]+=1
            else:
                if sex not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][sex]=1
                else:
                    cancer_cya_oa[ctype]['oa'][sex]+=1
        else:
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']={}
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']={}
            if '>' in age:
                if sex not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][sex]=1
                else:
                    cancer_cya_oa[ctype]['oa'][sex]+=1
            elif age == '<18':
                if sex not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][sex]=1
                else:
                    cancer_cya_oa[ctype]['cya'][sex]+=1
            elif int(age)<41:
                if sex not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][sex]=1
                else:
                    cancer_cya_oa[ctype]['cya'][sex]+=1
            else:
                if sex not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][sex]=1
                else:
                    cancer_cya_oa[ctype]['oa'][sex]+=1

    r.close()
    
    #print(cancer_cya_oa)
    w=open(outputfile, 'w')
    w.write('Type\tDisease\tCYA_Female_ratio\tOA_Female_ratio\tTotal\tCYA_Female\tCYA_Male\tOA_Female\tOA_Male\n')
    for ctype in cancer_cya_oa:
        if 'Female' not in cancer_cya_oa[ctype]['cya']:cancer_cya_oa[ctype]['cya']['Female']=0
        if 'Male' not in cancer_cya_oa[ctype]['cya']:cancer_cya_oa[ctype]['cya']['Male']=0
        if 'Female' not in cancer_cya_oa[ctype]['oa']:cancer_cya_oa[ctype]['oa']['Female']=0
        if 'Male' not in cancer_cya_oa[ctype]['oa']:cancer_cya_oa[ctype]['oa']['Male']=0
        total=cancer_cya_oa[ctype]['cya']['Female']+cancer_cya_oa[ctype]['cya']['Male']+cancer_cya_oa[ctype]['oa']['Female']+cancer_cya_oa[ctype]['oa']['Male']
        try:
            w.write(c_total_type[ctype]+'\t'+ctype+'\t'+str(cancer_cya_oa[ctype]['cya']['Female']/(cancer_cya_oa[ctype]['cya']['Male']+cancer_cya_oa[ctype]['cya']['Female']))+'\t'+str(cancer_cya_oa[ctype]['oa']['Female']/(cancer_cya_oa[ctype]['oa']['Male']+cancer_cya_oa[ctype]['oa']['Female']))+'\t'+str(total)+'\t'+str(cancer_cya_oa[ctype]['cya']['Female'])+'\t'+str(cancer_cya_oa[ctype]['cya']['Male'])+'\t'+str(cancer_cya_oa[ctype]['oa']['Female'])+'\t'+str(cancer_cya_oa[ctype]['oa']['Male'])+'\n')
        except:
            if cancer_cya_oa[ctype]['cya']['Male']==0:
                atmp='NA'
            else:
                atmp=str(cancer_cya_oa[ctype]['cya']['Female']/(cancer_cya_oa[ctype]['cya']['Male']+cancer_cya_oa[ctype]['cya']['Female']))
            if cancer_cya_oa[ctype]['oa']['Male']==0:
                btmp='NA'
            else:
                btmp=str(cancer_cya_oa[ctype]['oa']['Female']/(cancer_cya_oa[ctype]['oa']['Male']+cancer_cya_oa[ctype]['oa']['Female']))
            w.write(c_total_type[ctype]+'\t'+ctype+'\t'+atmp+'\t'+btmp+'\t'+str(total)+'\t'+str(cancer_cya_oa[ctype]['cya']['Female'])+'\t'+str(cancer_cya_oa[ctype]['cya']['Male'])+'\t'+str(cancer_cya_oa[ctype]['oa']['Female'])+'\t'+str(cancer_cya_oa[ctype]['oa']['Male'])+'\n')

            #print('Both Female and Male number == 0!')
            #print(ctype)
    w.close()

def get_cancer_metastasis(inforfile, cancer_type, outputfile, c_total_type):
    r=open(inforfile)
    cancer_cya_oa = {}
    r.readline()
    for line in r.readlines():
        line=line.strip()
        age=line.split('\t')[2]
        cancer=line.split('\t')[7]
        metastasis=line.split('\t')[4] # Primary or metastasis delete not collected and others
        if metastasis not in ['Metastasis','Primary']:
            metastasis='unknown'

        if cancer not in cancer_type:
            #print('error', cancer)
            continue
        ctype=cancer_type[cancer]
        if ctype not in cancer_cya_oa:
            cancer_cya_oa[ctype]={}
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']={}
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']={}
            if '>' in age:
                if metastasis not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['oa'][metastasis]+=1
            elif age == '<18':
                if metastasis not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['cya'][metastasis]+=1
            elif int(age)<41:
                if metastasis not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['cya'][metastasis]+=1
            else:
                if metastasis not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['oa'][metastasis]+=1
        else:
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']={}
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']={}
            if '>' in age:
                if metastasis not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['oa'][metastasis]+=1
            elif age == '<18':
                if metastasis not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['cya'][metastasis]+=1
            elif int(age)<41:
                if metastasis not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['cya'][metastasis]+=1
            else:
                if metastasis not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][metastasis]=1
                else:
                    cancer_cya_oa[ctype]['oa'][metastasis]+=1

    r.close()
    
    #print(cancer_cya_oa)
    w=open(outputfile, 'w')
    w.write('Type\tDisease\tCYA_Metastasis/Primary_ratio\tOA_Metastasis/Primary_ratio\tTotal\tCYA_Metastasis\tCYA_Primary\tOA_Metastasis\tOA_Primary\tCYA_Primary/Metastasis_ratio\tOA_Primary/Metastasis_tatio\n')
    for ctype in cancer_cya_oa:
        if 'Metastasis' not in cancer_cya_oa[ctype]['cya']:cancer_cya_oa[ctype]['cya']['Metastasis']=0
        if 'unknown' not in cancer_cya_oa[ctype]['cya']:cancer_cya_oa[ctype]['cya']['unknown']=0
        if 'Primary' not in cancer_cya_oa[ctype]['cya']:cancer_cya_oa[ctype]['cya']['Primary']=0
        if 'Metastasis' not in cancer_cya_oa[ctype]['oa']:cancer_cya_oa[ctype]['oa']['Metastasis']=0
        if 'unknown' not in cancer_cya_oa[ctype]['oa']:cancer_cya_oa[ctype]['oa']['unknown']=0
        if 'Primary' not in cancer_cya_oa[ctype]['oa']:cancer_cya_oa[ctype]['oa']['Primary']=0
        total=cancer_cya_oa[ctype]['cya']['Metastasis']+cancer_cya_oa[ctype]['cya']['Primary']+cancer_cya_oa[ctype]['oa']['Metastasis']+cancer_cya_oa[ctype]['oa']['Primary']+cancer_cya_oa[ctype]['cya']['unknown']+cancer_cya_oa[ctype]['oa']['unknown']
        try:
            w.write(c_total_type[ctype]+'\t'+ctype+'\t'+str(cancer_cya_oa[ctype]['cya']['Metastasis']/(cancer_cya_oa[ctype]['cya']['Primary']+cancer_cya_oa[ctype]['cya']['Metastasis']+cancer_cya_oa[ctype]['cya']['unknown']))+'\t'+str(cancer_cya_oa[ctype]['oa']['Metastasis']/(cancer_cya_oa[ctype]['oa']['Primary']+cancer_cya_oa[ctype]['oa']['Metastasis']+cancer_cya_oa[ctype]['oa']['unknown']))+'\t'+str(total)+\
                    '\t'+str(cancer_cya_oa[ctype]['cya']['Metastasis'])+'\t'+str(cancer_cya_oa[ctype]['cya']['Primary'])+\
                    '\t'+str(cancer_cya_oa[ctype]['oa']['Metastasis'])+'\t'+str(cancer_cya_oa[ctype]['oa']['Primary'])+'\t'+\
                    str(cancer_cya_oa[ctype]['cya']['Primary']/(cancer_cya_oa[ctype]['cya']['Primary']+cancer_cya_oa[ctype]['cya']['Metastasis']+cancer_cya_oa[ctype]['cya']['unknown']))+'\t'+str(cancer_cya_oa[ctype]['oa']['Primary']/(cancer_cya_oa[ctype]['oa']['Primary']+cancer_cya_oa[ctype]['oa']['Metastasis']+cancer_cya_oa[ctype]['oa']['unknown']))+'\n')
        except:
            w.write(c_total_type[ctype]+'\t'+ctype+'\tNA\tNA'+'\t'+str(total)+\
                    '\t'+str(cancer_cya_oa[ctype]['cya']['Metastasis'])+'\t'+str(cancer_cya_oa[ctype]['cya']['Primary'])+\
                    '\t'+str(cancer_cya_oa[ctype]['oa']['Metastasis'])+'\t'+str(cancer_cya_oa[ctype]['oa']['Primary'])+'\tNA\tNA\n')
    w.close()

def get_cancer_race(inforfile, cancer_type, outputfile, c_total_type):
    r=open(inforfile)
    cancer_cya_oa = {}
    races=[]
    r.readline()
    for line in r.readlines():
        line=line.strip()
        age=line.split('\t')[2]
        cancer=line.split('\t')[7]
        race1=line.split('\t')[10]  
        race2=line.split('\t')[11]
        race=race1+'_'+race2
        if 'Not' in race:continue
        if 'Other' in race:continue
        if 'Unknown' in race:continue
        if 'Pacific' in race:continue
        if race2=='Spanish/Hispanic':race=race2
        if race not in races:races.append(race)
        if cancer not in cancer_type:
            #print('error', cancer)
            continue
        ctype=cancer_type[cancer]
        if ctype not in cancer_cya_oa:
            cancer_cya_oa[ctype]={}
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']={}
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']={}
            if '>' in age:
                if race not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][race]=1
                else:
                    cancer_cya_oa[ctype]['oa'][race]+=1
            elif age == '<18':
                if race not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][race]=1
                else:
                    cancer_cya_oa[ctype]['cya'][race]+=1
            elif int(age)<41:
                if race not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][race]=1
                else:
                    cancer_cya_oa[ctype]['cya'][race]+=1
            else:
                if race not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][race]=1
                else:
                    cancer_cya_oa[ctype]['oa'][race]+=1
        else:
            if 'oa' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['oa']={}
            if 'cya' not in cancer_cya_oa[ctype]:cancer_cya_oa[ctype]['cya']={}
            if '>' in age:
                if race not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][race]=1
                else:
                    cancer_cya_oa[ctype]['oa'][race]+=1
            elif age == '<18':
                if race not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][race]=1
                else:
                    cancer_cya_oa[ctype]['cya'][race]+=1
            elif int(age)<41:
                if race not in cancer_cya_oa[ctype]['cya']:
                    cancer_cya_oa[ctype]['cya'][race]=1
                else:
                    cancer_cya_oa[ctype]['cya'][race]+=1
            else:
                if race not in cancer_cya_oa[ctype]['oa']:
                    cancer_cya_oa[ctype]['oa'][race]=1
                else:
                    cancer_cya_oa[ctype]['oa'][race]+=1

    r.close()
    
    #print(cancer_cya_oa)
    w=open(outputfile, 'w')
    w.write('Type\tDisease')
    for race in races:
        w.write('\tCYA_'+race+'\tOA_'+race)
    w.write('\n')

    for ctype in cancer_cya_oa:
        w.write(c_total_type[ctype]+'\t'+ctype+'\t')
        for race in races:
            if race not in cancer_cya_oa[ctype]['cya']:
                w.write('0\t')
            else:
                w.write(str(cancer_cya_oa[ctype]['cya'][race])+'\t')
            if race not in cancer_cya_oa[ctype]['oa']:
                w.write('0\t')
            else:
                w.write(str(cancer_cya_oa[ctype]['oa'][race])+'\t')
        w.write('\n')
    w.close()

def get_cancer_type(cancertypefile):
    r=open(cancertypefile)
    c_total={}
    for line in r.readlines():
        line=line.strip()
        c_total[line.split('\t')[1]]=line.split('\t')[0]
    return(c_total)

def main():
    dividefile='hematologic_divided.txt'
    inforfile='patient_information.txt'
    outputfile='results/result_cancer_age.txt'
    sex_outputfile='results/result_cancer_sex.txt'
    metastasis_outputfile='results/result_cancer_metastasis.txt'
    cancertypefile='results/cancer_type.txt'
    race_outputfile='results/result_cancer_race.txt'
    
    cancer_type=get_cancer_infor(dividefile)
    c_total_type=get_cancer_type(cancertypefile)
    get_cancer_age(inforfile, cancer_type, outputfile, c_total_type)
    get_cancer_sex(inforfile, cancer_type, sex_outputfile, c_total_type)
    get_cancer_metastasis(inforfile, cancer_type, metastasis_outputfile, c_total_type)
    get_cancer_race(inforfile, cancer_type, race_outputfile, c_total_type)

if __name__=='__main__':
    main()

        
    
