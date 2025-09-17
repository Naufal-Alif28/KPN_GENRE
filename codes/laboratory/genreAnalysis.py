#!/usr/bin/env python3

#import
from cobra.io import read_sbml_model
import riptide as rp
import csv


def RichMedium(model):
    #components
    LB_5 = {'EX_glc__D_e', 'EX_ala_B_e', 'EX_ala__D_e', 'EX_ala__L_e', 'EX_arg__L_e', 'EX_asp__L_e', 'EX_cys__L_e',
        'EX_glu__L_e', 'EX_gly_e', 'EX_his__L_e', 'EX_ile__L_e', 'EX_leu__L_e', 'EX_lys__L_e', 'EX_met__L_e',
        'EX_pro__L_e', 'EX_thr__L_e', 'EX_tyr__L_e', 'EX_phe__L_e', 'EX_ser__D_e', 'EX_ser__L_e', 'EX_trp__L_e',
        'EX_val__L_e', 'EX_pnto__R_e', 'EX_nac_e', 'EX_ins_e', 'EX_hxan_e', 'EX_dcyt_e', 'EX_thymd_e',
        'EX_ura_e', 'EX_uri_e', 'EX_dad_2_e', 'EX_adn_e', 'EX_fru_e', 'EX_gal_e'}
    LB_1000 = {'EX_na1_e', 'EX_cl_e', 'EX_so4_e', 'EX_k_e', 'EX_pi_e', 'EX_ca2_e', 'EX_mg2_e', 'EX_zn2_e', 'EX_aso3_e',
            'EX_cd2_e', 'EX_hg2_e', 'EX_co2_e', 'EX_cobalt2_e', 'EX_cu2_e', 'EX_fe2_e', 'EX_fe3_e', 'EX_mn2_e',
            'EX_mobd_e', 'EX_tungs_e', 'EX_ni2_e'}
    LB_100 = {'EX_h_e', 'EX_h2o_e'}
    LB_18_5 = {'EX_o2_e'}
    LB_0_01 = {'EX_cbl1_e'}
    #reset boundary for all components
    for x in model.boundary:
        x.lower_bound = 0.0
        x.upper_bound = 0.0
    #boundary map (LB_,lower_boundary,upper_boundary)
    boundary_map = [
        ("5",-5.0,1000.0),("1000",-1000.0,1000.0),("100",-100.0,1000.0),
        ("18_5",-18.5,1000.0),("0_01",-0.01,1000.0)
    ]
    #changing boundary
    for rxn in model.boundary:
        for x in boundary_map:
            if rxn.id in locals()[f"LB_{x[0]}"]:
                rxn.lower_bound = x[1]
                rxn.upper_bound = x[2]
    
    return model

def main():
    #start
    print("start genreAnalysis.py")

    #SRA list maker
    sra_list = []
    with open("sra_list.txt",mode='r') as file:
        for sra in file:
            sra = sra.rstrip()
            sra_list.append(sra)

    with open("expressionAnalysis.csv",mode='r') as file:
        exp_csvobj = csv.reader(file, delimiter=',')
        lines = list(exp_csvobj)
    for sra in sra_list:
        gene_id_list = []
        for line in lines[1:]:
            gene_id_list.append(line[0])
        sra_idx = lines[0].index(f'{sra}')
        exp_list = []
        for line in lines[1:]:
            exp_list.append(line[1]) 
        rows = [['gene_id','transcript_abund']]
        for i in range(len(gene_id_list)):
            rows.append([gene_id_list[i],exp_list[i]])
        with open(f"genreAnalysis/transcriptome_{sra}.tsv",mode='w',newline='') as file:
            tsv_out = csv.writer(file,delimiter='\t')
            tsv_out.writerows(rows)

    #MODEL CONTEXTUALISATION
    #read reference model
    temp_model = read_sbml_model("../iYL1228.xml")
    #growth medium definition
    model = RichMedium(temp_model)
    #contextualisation
    for sra in sra_list:
        trans_dict = rp.read_transcription_file(f"genreAnalysis/transcriptome_{sra}.tsv", header=True)
        model_context = rp.contextualize(model=model, transcriptome=trans_dict, fraction=0.75, threshold=1e-6, silent=True) 
        rp.save_output(riptide_obj=model_context, path=f"genreAnalysis/{sra}", file_type="sbml")

    #end
    print("end genreAnalysis.py")

if __name__ == '__main__':
    main()