import sys, os
#os.chdir('C:/Users/Li/Desktop/CDAP/scripts/')
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

map_file=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Testing\\4_CPTAC2_Ovarian_Prospective_Collection_PNNL\\Supporting_Documents\\Proteome_Mapping_PNNL.txt')
PNNL_scheme = dict(zip(map_file['Specimen Label'], map_file['Participant ID']))

prospective_proteome_summary_report=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Testing\\4_CPTAC2_Ovarian_Prospective_Collection_PNNL\\SummaryReports\\CPTAC2_Ovarian_Prospective_Collection_PNNL_Proteome.tmt10.tsv')
#table=prospective_proteome_summary_report

def save_reorder_column(path,table,name):
    table=table.groupby(['Gene_ID'],as_index=False).mean()
    table.sort_values('Gene_ID', inplace=True) #sort by Gene_ID
    table.sort_index(axis=1, inplace=True) #sort column names
    cols = list(table)
    # move the column to head of list using index, pop and insert
    cols.insert(0, cols.pop(cols.index('Gene_ID')))
    table = table.ix[:, cols]
    table.to_csv(path+name+'.txt',sep='\t',index=False)
    
def clean_proteome(table,name):
    table = pd.concat([table[[col for col in table.columns if 'Unshared' in col]], table['Gene']], axis=1)#keep only Unshared Version of the quantitation
    ##remove 19 normal samples 
    table=table.ix[3:] #remove first three rows, mean, median, std
    table_new_colnames=[]
    for i in list(table.columns.values):
        table_new_colnames.append(i.split(' ')[0])
    table.columns = table_new_colnames
    table = table.rename(columns=PNNL_scheme)    
    cols = [c for c in table.columns if (c[:3] != 'JHU') and (c[:1]!='N')] #remove JHU.QC samples and normal samples
    table=table[cols] #84 PNNL cancer proteome
    table=table.rename(columns={'Gene': 'Gene_ID'})
    
    table_filtered=table.dropna(thresh=5) #remove rows with no more than 5 no-nan values
    save_reorder_column('C:/Users/Li/Desktop/',table,name+'_all_gene')
    save_reorder_column('C:/Users/Li/Desktop/',table_filtered,name+'_filtered')

clean_proteome(prospective_proteome_summary_report,'prospective_proteome_PNNL')

prospecitve_phospho_summary_report=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Testing\\4_CPTAC2_Ovarian_Prospective_Collection_PNNL\\SummaryReports\\CPTAC2_Ovarian_Prospective_Collection_PNNL_Phosphoproteome.phosphosite.tmt10.tsv')

def clean_phospho(table,name):
    table.drop(['Peptide','Organism'],inplace=True,axis=1) #remove Peptide, Organism columns
    try:
        table.drop(['263d3f-I Log Ratio','blcdb9-I Log Ratio','c4155b-C Log Ratio'],inplace=True,axis=1) #drop three normal samples in breast cancer data  
    except ValueError:
        pass  # do nothing!    
    #table.drop(table.columns[[1,2,3]], axis=1, inplace=True) #remove 3 normal sample
    table=table.ix[3:] #remove first three rows, mean, median, std
    table_new_colnames=[]
    for i in list(table.columns.values):
        table_new_colnames.append(i.split(' ')[0])
    table.columns = table_new_colnames
    table = table.rename(columns=PNNL_scheme)    
    cols = [c for c in table.columns if (c[:3] != 'JHU') and (c[:1]!='N')] #remove JHU.QC samples and normal samples
    table=table[cols] #84 PNNL cancer proteome

    table['Gene_ID'] = table[['Gene', 'Phosphosite']].apply(lambda x: '.'.join(x), axis=1)
    table.drop(['Gene','Phosphosite'],inplace=True,axis=1) #remove Peptide, Organism columns; 
    
    table_filtered=table.dropna(thresh=5) #remove rows with no more than 5 no-nan values 
    save_reorder_column('C:/Users/Li/Desktop/',table,name+'_all_gene')
    save_reorder_column('C:/Users/Li/Desktop/',table_filtered,name+'_filtered')

clean_phospho(prospecitve_phospho_summary_report,'prospective_phospho_PNNL')


##CDAP summary table:
#full size of PNNL phospho 45622*83, filtered: 45284*， 5QCs, 20 Normal Tissue
#full size of PNNL proteome 8818*83, filtered:8818*83 , 5QC, 20 Normal Tissue
#7 sample with proteome measurement missing in phospho: {'02OV029', '02OV044', '04OV001', '04OV004', '04OV008', '04OV049', '17OV001'}
#20 normal tissues 


##read in prospective filtered proteome/phospho data processed by CDAP
pros_phos_raw=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Testing\\4_CPTAC2_Ovarian_Prospective_Collection_PNNL\\prospective_phospho_PNNL_filtered.txt')
pros_proteome_raw=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Testing\\4_CPTAC2_Ovarian_Prospective_Collection_PNNL\\prospective_proteome_PNNL_filtered.txt')

##read in prospective RNA-seq and CNA data processed by GDC, reanalyzed by xuya and zhi
pros_rna_raw=pd.read_table('C:\\Dream_data\\data_for_synapse\\testing\\sub2\\complete_data\\RNA_seq_ova_prospective_fpkm_xuya.txt')
pros_cna_raw=pd.read_table('C:\\Dream_data\\data_for_synapse\\testing\\sub2\\complete_data\\CNA_WXS_ovarian_prospective_median.txt')

pros_rna_filtered=pros_rna_raw.loc[np.sum(pros_rna_raw.iloc[:,]>0, axis=1)[lambda x: x>=5].index] #filter samples with less than 5 non-zero values

len(set(pros_rna_raw['Gene_ID'])-set(pros_rna_filtered['Gene_ID']))

pros_cna_filtered=pros_cna_raw.loc[np.sum(pros_cna_raw.iloc[:,]>0, axis=1)[lambda x: x>=5].index] #filter samples with less than 5 non-zero values

##find common samples in all platform (82 samples) RNA-seq of 01OV029 is disqualified by GDC)
PNNL_ova_shared_names=sorted(list(set(pros_proteome_raw.columns.values) & set(pros_phos_raw.columns.values)&set(pros_rna_raw.columns.values)&set(pros_cna_raw.columns.values)))
PNNL_ova_shared_names_proteins=sorted(list(set(pros_proteome_raw.columns.values) & set(pros_phos_raw.columns.values)))

pros_proteome_common_sample=pros_proteome_raw[PNNL_ova_shared_names]
pros_phos_common_sample=pros_phos_raw[PNNL_ova_shared_names]
#pros_rna_common_sample=pros_rna_raw[PNNL_ova_shared_names]
#pros_cna_common_sample=pros_cna_raw[PNNL_ova_shared_names]
pros_cna_common_sample=pros_cna_filtered[PNNL_ova_shared_names]
pros_rna_common_sample=pros_rna_filtered[PNNL_ova_shared_names]

##read in retro proteome/phospho data and RNA-seq/Microarray, CNA data
retro_ova_CNA_raw=pd.read_table('C:/DARPA/CPTAC_Ovarian/Training_data(CDAP)/cna_ovarian',sep=',')
retro_ova_CNA=retro_ova_CNA_raw.rename(columns={'Unnamed: 0': 'Gene_ID'})

retro_ova_rna=pd.read_table('C:/DARPA/CPTAC_Ovarian/Training_data(CDAP)/rna_ovarian',sep=',')
retro_ova_rna=retro_ova_rna.rename(columns={'Unnamed: 0': 'Gene_ID'})

retro_ova_rna_seq=pd.read_table('C:/DARPA/CPTAC_Ovarian/OV.uncv2.mRNAseq_RSEM_Z_Score.txt',sep=',')
retro_ova_rna_seq=retro_ova_rna_seq.rename(columns={'Gene': 'Gene_ID'})
retro_ova_rna_seq=retro_ova_rna_seq.drop(retro_ova_rna_seq.columns[[0]], axis=1)  
retro_ova_rna_seq=retro_ova_rna_seq.groupby(['Gene_ID'],as_index=False).mean()

retro_ova_pro_JHU_filtered=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Training\\cleaned_data\\retrospective_ova_proteome_JHU_filtered.txt')
retro_ova_pro_PNNL_filtered=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Training\\cleaned_data\\retrospective_ova_proteome_PNNL_filtered.txt')
retro_ova_phospho_PNNL_filtered=pd.read_table('C:\\Dream_data\\CDAP\\CDAP_data(Paul_Nathan_version)\\Training\\cleaned_data\\retrospective_ova_phospho_filtered.txt')
##read in retro proteome/phospho data and RNA-seq/Microarray, CNA data

##find common gene_ids shared by retro and prospective 

def common_gene_names(*args):
    result = set(args[0])
    for s in args[1:]:
        result.intersection_update(s)
    return sorted(list(result))

def save_common_gene(path,table,commom_gene_list,name,type):
    table=table.groupby(['Gene_ID'],as_index=False).mean()
    table_commongene=table[table['Gene_ID'].isin(commom_gene_list)]
    table_commongene.sort_values('Gene_ID', inplace=True)
    cols = list(table_commongene)
    cols.insert(0, cols.pop(cols.index('Gene_ID')))
    table_commongene = table_commongene.ix[:, cols]
    if type=='prospective':
        table_commongene.columns = ['Participant_'+str(i-1) if 2 <= i  else x for i, x in enumerate(table_commongene.columns, 1)]
    table_commongene.to_csv(path+name+'_common_gene_'+str(len(commom_gene_list))+'.txt',sep='\t',index=False)

ova_cna_common_gene_names=common_gene_names(pros_cna_common_sample['Gene_ID'],retro_ova_CNA['Gene_ID']) #11859
ova_rna_common_gene_names=common_gene_names(pros_rna_common_sample['Gene_ID'],retro_ova_rna['Gene_ID'],retro_ova_rna_seq['Gene_ID']) #15121
ova_pro_common_gene_names=common_gene_names(pros_proteome_common_sample['Gene_ID'],retro_ova_pro_PNNL_filtered['Gene_ID'],retro_ova_pro_JHU_filtered['Gene_ID']) #7229 gene passed nonmissing rate shared by both PNNL retro and prospective data (7061 shared by retro PNNL JUH and PNNL prospective)
ova_phospho_common_gene_names=common_gene_names(pros_phos_common_sample['Gene_ID'],retro_ova_phospho_PNNL_filtered['Gene_ID']) #12237 phosphosites

save_common_gene('C:/Users/Li/Desktop/',pros_cna_common_sample, ova_cna_common_gene_names,'prospective_ova_CNA_sort','prospective')
save_common_gene('C:/Users/Li/Desktop/',pros_rna_common_sample, ova_rna_common_gene_names,'prospective_ova_rnaseq_sort','prospective')

save_common_gene('C:/Users/Li/Desktop/',pros_proteome_common_sample, ova_pro_common_gene_names,'prospective_ova_proteome_sort','prospective') #7229 gene passed nonmissing rate shared by both PNNL retro and prospective data
save_common_gene('C:/Users/Li/Desktop/',pros_phos_common_sample, ova_phospho_common_gene_names,'prospective_ova_phospho_sort','prospective')

save_common_gene('C:/Users/Li/Desktop/',retro_ova_CNA, ova_cna_common_gene_names,'retrospective_ova_CNA_sort','retro')
save_common_gene('C:/Users/Li/Desktop/',retro_ova_rna, ova_rna_common_gene_names,'retrospective_ova_array_sort','retro')
save_common_gene('C:/Users/Li/Desktop/',retro_ova_rna_seq, ova_rna_common_gene_names,'retrospective_ova_rna_seq_sort','retro')

save_common_gene('C:/Users/Li/Desktop/',retro_ova_pro_JHU_filtered, ova_pro_common_gene_names,'retrospective_ova_JHU_proteome_sort','retro') #7229 gene passed nonmissing rate shared by both PNNL retro and prospective data
save_common_gene('C:/Users/Li/Desktop/',retro_ova_pro_PNNL_filtered, ova_pro_common_gene_names,'retrospective_ova_PNNL_proteome_sort','retro') #7229 gene passed nonmissing rate shared by both PNNL retro and prospective data
save_common_gene('C:/Users/Li/Desktop/',retro_ova_phospho_PNNL_filtered, ova_phospho_common_gene_names,'retrospective_ova_phospho_sort','retro')

sum(pros_rna_filtered.std(axis=1)!=0)
#output processed data. update Mi's summary report
#send out group email to everyone