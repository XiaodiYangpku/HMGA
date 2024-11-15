#！/usr/bin/env python
#-*-coding:utf-8 -*-

import os,re
### select hematologic tumor patients from data_clinical_patients.txt and combine data_clinical_patient.txt to obtain DOB/DOD ###
### the 7 coloumn ###
# Blastic Plasmacytoid Dendritic Cell Neoplasm
# Blood Cancer
# Blood Cancer, NOS
# B-Lymphoblastic Leukemia/Lymphoma
# Hodgkin Lymphoma
# Leukemia
# Lymphatic Cancer, NOS
# Mastocytosis
# Mature B-Cell Neoplasms
# Mature T and NK Neoplasms
# Myelodysplastic/Myeloproliferative Neoplasms
# Myelodysplastic Syndromes
# Myeloid Neoplasms with Germ Line Predisposition
# Myeloproliferative Neoplasms
# Non-Hodgkin Lymphoma
# Posttransplant Lymphoproliferative Disorders
# T-Lymphoblastic Leukemia/Lymphoma

def get_patientid_status(patientfile):
    patientid_status={}
    with open(patientfile) as r:
        for line in r.readlines():
            line=line.strip('\n')
            if line[0]=='#':continue
            patientid=line.split('\t')[0]
            sex=line.split('\t')[1]
            race=line.split('\t')[2]
            dob_last_contact=line.split('\t')[5]
            dob_dod=line.split('\t')[6]
            status=line.split('\t')[8]
            if dob_last_contact=='Unknown':continue
            if dob_dod=='Unknown':continue
            if 'Not Collected' in dob_last_contact and '<' in dob_dod:continue
            if '<' in dob_last_contact and 'Not Collected' in dob_dod:continue
            if 'Not' in dob_last_contact and 'Not' in dob_dod:continue
            if dob_last_contact=='' and dob_dod=='':continue
            patientid_status[patientid]='\t'.join(line.split('\t')[1:])
    return(patientid_status)

def get_patientid_cancer(samplefile):
    cancers=['CANCER_TYPE','Blastic Plasmacytoid Dendritic Cell Neoplasm','Blood Cancer','Blood Cancer, NOS','B-Lymphoblastic Leukemia/Lymphoma','Hodgkin Lymphoma','Leukemia','Lymphatic Cancer, NOS','Mastocytosis','Mature B-Cell Neoplasms','Mature T and NK Neoplasms','Myelodysplastic/Myeloproliferative Neoplasms','Myelodysplastic Syndromes','Myeloid Neoplasms with Germ Line Predisposition','Myeloproliferative Neoplasms','Non-Hodgkin Lymphoma','Posttransplant Lymphoproliferative Disorders','T-Lymphoblastic Leukemia/Lymphoma']
    patientid_cancer={}
    with open(samplefile) as r:
        for line in r.readlines():
            line=line.strip('\n')
            if line[0]=='#':continue
            patientid=line.split('\t')[0]
            sampleid=line.split('\t')[1]
            age=line.split('\t')[2]
            if age=='Unknown':continue
            oncotree_code=line.split('\t')[3]
            sample_type=line.split('\t')[4]
            cancer_type=line.split('\t')[6]
            detailed_cancer=line.split('\t')[7]
            detailed_sample=line.split('\t')[8]
            if cancer_type not in cancers:continue
            if patientid not in patientid_cancer:
                patientid_cancer[patientid]='\t'.join(line.split('\t')[1:])
            else:
                pass
                #print(line)
    return(patientid_cancer)

def get_patientid_infor(patientid_status,patientid_cancer,outputfile):
    w=open(outputfile,'w')
    for patientid in patientid_cancer:
        if patientid not in patientid_status:continue
        w.write(patientid+'\t'+patientid_cancer[patientid]+'\t'+patientid_status[patientid]+'\n')
    w.close()

def main():
    samplefile='genie-data/data_clinical_sample.txt'
    patientfile='genie-data/data_clinical_patient.txt'
    outputfile='patient_information.txt'

    patientid_status=get_patientid_status(patientfile)
    patientid_cancer=get_patientid_cancer(samplefile)
    get_patientid_infor(patientid_status,patientid_cancer,outputfile)

if __name__=='__main__':
    main()




