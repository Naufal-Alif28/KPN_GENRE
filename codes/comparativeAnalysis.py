#!/usr/bin/env python3

#start
print("start comparativeAnalysis.py")

#IMPORT
from cobra.io import read_sbml_model
from cobra.flux_analysis.variability import find_essential_reactions
import pandas as pd
#from concurrent import futures

#FUNCTIONS
#BinSet: 
#function to create a set of reactions (isolate_rxn_set), reactions existance binary dict (isolate_bin_dict)
#a list of number of reactions per isolate group (rxn_num), and a list of number of reactions per SRA (sra_rxn)
def BinSet(isolate:str, ref_rxn_list):
    sra_list = []
    with open(f"{isolate}/sra_list.txt", mode='r') as file:
        for sra in file:
            sra = sra.rstrip()
            sra_list.append(sra)

    #create inital dictionary for storing rxn existance bin
    rxn_bin_dict = {}

    isolate_rxn_set = set()
    sra_rxn = []
    rxn_num = [0 for i in range(len(ref_rxn_list))]
    for sra in sra_list:
        model = read_sbml_model(f"{isolate}/genreAnalysis/{sra}/model.sbml")
        sra_rxn_set = set()
        #adding reactions to the set (excluding duplicate)
        for rxn in [rxn.id for rxn in model.reactions]:
            sra_rxn_set.add(rxn)
        isolate_rxn_set = isolate_rxn_set | sra_rxn_set #union --> return
        #adding rxn existance bin to rxn_dict
        rxn_bin_list = []
        idx = 0
        for rxn in ref_rxn_list:
            if rxn in isolate_rxn_set:
                rxn_bin_list.append(1)
                rxn_num[idx] += 1 #--> return
            else:
                rxn_bin_list.append(0)
            idx += 1
        rxn_bin_dict[f'{sra}'] = rxn_bin_list #--> return
        sra_rxn.append(len(model.reactions)) #--> return

    return isolate_rxn_set, rxn_bin_dict, rxn_num, sra_rxn

#ExcelOutput:
#Excel output function
def ExcelOutput(prefix:str, bin_dict:dict, rxn_num:dict, sra_rxn:dict, is_new=bool):
    if is_new:
        init_mode = 'w'
    elif not is_new:
        init_mode = 'a'
    with pd.ExcelWriter("comparativeAnalysisOutput.xlsx",mode=init_mode) as comparativeWriter:
        #reactions binary state csv output (merge_bin_dict)
        df = pd.DataFrame(data=bin_dict)
        df.to_excel(comparativeWriter,sheet_name=f"{prefix}_bin")
    with pd.ExcelWriter("comparativeAnalysisOutput.xlsx",mode='a') as comparativeWriter:
        #output rxn_num_dict
        df = pd.DataFrame(data=rxn_num)
        df.to_excel(comparativeWriter,sheet_name=f"{prefix}_rxn_num")
        # number of reactions per SRA (sra_rxn_dict)
        df = pd.DataFrame(data=sra_rxn)
        df.to_excel(comparativeWriter,sheet_name=f"{prefix}_sra_rxn")

def main():
    #PROGRAM
    #CONTEXTUALISED MODELS ANALYSIS
    #define a list of SRAs (sra_list) and sets of clinical and laboratory SRAs
    sra_list = []
    sra_clinical_set = set()
    sra_laboratory_set = set()
    with open("clinical/sra_list.txt", mode='r') as file:
        for sra in file:
            sra = sra.rstrip()
            sra_list.append(sra)
            sra_clinical_set.add(sra)
    with open("laboratory/sra_list.txt", mode='r') as file:
        for sra in file:
            sra = sra.rstrip()
            sra_list.append(sra)
            sra_laboratory_set.add(sra)
    #read reference model
    model = read_sbml_model('iYL1228.xml')
    #create list of reference reactions
    ref_rxns = [rxn.id for rxn in model.reactions] 
    #create reactions existance binary dictionary
    clinical_rxn_set, clinical_bin_dict, clinical_rxn_num, clinical_sra_rxn = BinSet(isolate='clinical',ref_rxn_list=ref_rxns)
    laboratory_rxn_set, laboratory_bin_dict, laboratory_rxn_num, laboratory_sra_rxn = BinSet(isolate='laboratory',ref_rxn_list=ref_rxns)
    #merging and indexing dictionaries
    merge_bin_dict = {
        'rxn_id' : ref_rxns
    } | clinical_bin_dict | laboratory_bin_dict

    rxn_num_dict = {
            'rxn_id' : [rxn.id for rxn in model.reactions],
            'clinical_rxn_num' : clinical_rxn_num,
            'laboratory_rxn_num' : laboratory_rxn_num
        }
    sra_rxn_dict = {
        'SRA' : sra_list,
        'rxn_num' : clinical_sra_rxn + laboratory_sra_rxn
    }
    #Excel output for contextualised models analysis ('all')
    ExcelOutput(prefix="all", bin_dict=merge_bin_dict, rxn_num=rxn_num_dict, sra_rxn=sra_rxn_dict, is_new=True)

    #REACTIONS ESSENTIALITY ANALYSIS
    #biomass objective value optimalisation of reference model
    optimum_growth = model.optimize().objective_value

    #update merge_bin_dict and rxn_set
    for isolate in ["clinical","laboratory"]:
        for sra in list(locals()[f'sra_{isolate}_set']):
            model = read_sbml_model(f"{isolate}/genreAnalysis/{sra}/model.sbml")
            merge_bin_dict[f'{sra}'] = [0 for i in range(len(merge_bin_dict['rxn_id']))]
            essential_rxn_set = find_essential_reactions(model)
            essential_rxn_list = [rxn.id for rxn in essential_rxn_set]
            for rxn in essential_rxn_list:
                idx = merge_bin_dict['rxn_id'].index(rxn)
                merge_bin_dict[f'{sra}'][idx] = 1

    #update sra_rxn 
    for sra in sra_list:
        idx = sra_rxn_dict['SRA'].index(sra)
        sra_rxn_dict['rxn_num'][idx] = sum(merge_bin_dict[f'{sra}'])
    #update rxn_num 
    for rxn in rxn_num_dict['rxn_id']:
        idx1 = merge_bin_dict['rxn_id'].index(rxn)
        idx2 = rxn_num_dict['rxn_id'].index(rxn)
        rxn_num_dict['clinical_rxn_num'][idx2] = 0
        rxn_num_dict['laboratory_rxn_num'][idx2] = 0
        for sra in sra_list:
            if sra in sra_clinical_set:
                rxn_num_dict['clinical_rxn_num'][idx2] += merge_bin_dict[f'{sra}'][idx1]
            elif sra in sra_laboratory_set:
                rxn_num_dict['laboratory_rxn_num'][idx2] += merge_bin_dict[f'{sra}'][idx1]
    #update clinical & laboratory rxn set
    for i in range(len(rxn_num_dict['rxn_id'])):
        rxn = rxn_num_dict['rxn_id'][i]
        if (rxn in clinical_rxn_set) and (rxn_num_dict['clinical_rxn_num'][i] == 0):
            clinical_rxn_set.remove(rxn)
        elif (rxn in laboratory_rxn_set) and (rxn_num_dict['laboratory_rxn_num'][i] == 0):
            laboratory_rxn_set.remove(rxn)
    #Excel output for pruned models with essential reactions only ('essential')
    ExcelOutput(prefix="essential", bin_dict=merge_bin_dict, rxn_num=rxn_num_dict, sra_rxn=sra_rxn_dict, is_new=False)

    #CORE REACTIONS SUBTRACTION
    #leaving reactions unique to each isolate
    #(I = C ∩ L)
    intersect_rxn_set = clinical_rxn_set.intersection(laboratory_rxn_set) #(I = C ∩ L)
    #update bin_dict
    pop_idx_list = []
    for rxn in merge_bin_dict['rxn_id']:
        if rxn in intersect_rxn_set:
            idx = merge_bin_dict['rxn_id'].index(rxn)
            pop_idx_list.append(idx)
    for i in sorted(pop_idx_list, reverse=True):
        merge_bin_dict['rxn_id'].pop(i)
    for sra in sra_list:
        for i in sorted(pop_idx_list, reverse=True):
            merge_bin_dict[f'{sra}'].pop(i)
    #update rxn_num
    for i in sorted(pop_idx_list, reverse=True):
        rxn_num_dict['rxn_id'].pop(i)
        rxn_num_dict['clinical_rxn_num'].pop(i)
        rxn_num_dict['laboratory_rxn_num'].pop(i)
    #update sra_rxn 
    for sra in sra_list:
        idx = sra_rxn_dict['SRA'].index(sra)
        sra_rxn_dict['rxn_num'][idx] = sum(merge_bin_dict[f'{sra}'])
    #Excel output updated for models with unique/accessory reactions ('unique')
    ExcelOutput(prefix="unique", bin_dict=merge_bin_dict, rxn_num=rxn_num_dict, sra_rxn=sra_rxn_dict, is_new=False)

    #COMMON REACTION IDENTIFICATION (55% share threshold)
    ids = rxn_num_dict['rxn_id']
    N = len(sra_clinical_set)
    common_rxn_clinical = []
    for i in range(len(ids)):
        if rxn_num_dict['clinical_rxn_num'][i] >= 0.55*N:
            common_rxn_clinical.append(ids[i])
    N = len(sra_laboratory_set)
    common_rxn_laboratory = []
    for i in range(len(ids)):
        if rxn_num_dict['laboratory_rxn_num'][i] >= 0.55*N:
            common_rxn_laboratory.append(ids[i])
    #output for summary in Excel
    with pd.ExcelWriter("comparativeAnalysisOutput.xlsx",mode='a') as comparativeWriter:
        df = pd.DataFrame(data={'clinical_common_rxn' : common_rxn_clinical})
        df.to_excel(comparativeWriter,sheet_name="clinical_common_rxn")
        df = pd.DataFrame(data={'laboratory_common_rxn' : common_rxn_laboratory})
        df.to_excel(comparativeWriter,sheet_name="laboratory_common_rxn")

    #end
    print("end comparativeAnalysis.py")

if __name__ == '__main__':
    main()