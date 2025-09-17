#!/usr/bin/env python3

def main():
    #start
    print("start expressionAnalysis")

    #import
    import csv
    import pandas as pd

    #read the list of SRA
    sra_list = []
    with open("sra_list.txt") as file:
        for sra in file:
            sra = sra.rstrip()
            sra_list.append(sra)

    #build initial dictionary
    num_genes = 4857 #number of genes in the sequence annotation refrence (.GTF file)
    gene_ids = []
    for i in range(1,num_genes+1):
        zeros = "0"*(5-len(str(i)))
        id = "KPN_"+zeros+str(i)
        gene_ids.append(id)
    exp_dict = {
        'gene_id' : gene_ids
    }

    #EXTRACTION AND NORMALIZATION
    for sra in sra_list:
        print(f"{sra}")
        with open(f"reads_annotate/{sra}_annotate.tsv",mode='r',newline='') as file:
            read_csvobj = csv.reader(file,delimiter='\t')
            lines = list(read_csvobj)
            #read and store number of reads (num_read) of each SRA
            num_read = [0 for i in range(num_genes)]
            for line in lines:
                try:
                    id_parse = line[0].split('_')
                    id = int(id_parse[1])
                    num_read[id-1] = int(line[6])
                except:
                    pass
            #read and store length of reads (len_read) of each SRA
            len_read = [0 for i in range(num_genes)]
            for line in lines:
                try:
                    id_parse = line[0].split('_')
                    id = int(id_parse[1])
                    len_read[id-1] = int(line[5])
                except:
                    pass
            #calculate normalised expression (norm_read) of each SRA (using TPM normalisation)
            norm_read = [0 for i in range(num_genes)]
            for i in range(num_genes):
                if len_read[i] != 0:
                    norm_read[i] = num_read[i]/(len_read[i]) #Divide the read counts by the length of each gene in kilobases (RPK, reads per kilobase )
            sf = sum(norm_read)/1000 #Count up all the RPK values in a sample and divide this number by 1,000,000 (sf, scaling factor)
            for i in range(num_genes):
                norm_read[i] = norm_read[i]/sf #Divide the RPK values by sf (TPM, Transcripts Per Kilobase Million)
        
        exp_dict[f'{sra}'] = norm_read

    #convert dictionary to Pandas DataFrame
    df = pd.DataFrame(data=exp_dict)
    #calculate mean and stdev of each gene's num_read
    num_col_core = len(df.columns)
    df['average'] = df.iloc[:,1:num_col_core].mean(axis=1)
    df['stdev'] = df.iloc[:,1:num_col_core].std(axis=1)
    #convert df to csv
    df.to_csv("expressionAnalysis.csv",index=False)

    #end
    print("end expressionAnalysis")

if __name__ == '__main__':
    main()