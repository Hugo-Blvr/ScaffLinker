import pandas as pd
import os

pd.set_option('display.max_rows', None)


def preprocess(paf_dir, NbMatch, IdSeq):

    paf_files = []
    for paf in os.listdir(paf_dir):paf_files.append(f"{paf_dir}/{paf}")
    df=pd.DataFrame()
    for paf_file in paf_files:    
        df_paf = pd.read_csv(paf_file,sep='\t',header=None).iloc[:, 0:12] # récupération des informations
        df_paf = df_paf[df_paf.iloc[:, 11] > 40].iloc[:, 0:11] # filtre sur la qualité de mapping
        paf_name = paf_file.split('/')[-1].split('.')[0].split('_')[-1] # détermination du nom de l'échantillon
        df_paf[0] = f"{paf_name}$" + df_paf[0]
        df = pd.concat([df,df_paf])
    
    df['IdSeq'] = df[9] / df[10]
    df.columns = ["Qname","Qlen","Qstart","Qstop","Strand","Tname","Tlen","Tstart","Tstop","NbMatch","NbBase","IdSeq"]
    df = df[["Qname","Qlen","Qstart","Qstop","Strand","Tname","Tlen","Tstart","Tstop","NbMatch","NbBase","IdSeq"]]
    
    df_filtre = df[(df["NbMatch"] >= NbMatch) & (df["IdSeq"] >= IdSeq)]
    df = df[(df["NbMatch"] < NbMatch) | (df["IdSeq"] < IdSeq)]
    
    return df.reset_index(drop=True), df_filtre.reset_index(drop=True)



def Run(paf_dir, NbMatch, IdSeq):  
    df, df_filtre=preprocess(paf_dir, NbMatch, IdSeq)

    print(df)

if __name__ == "__main__":
    paf_dir = "02_masked_paf_files/02_paf_files_Gd45"
    NbMatch = 5000
    IdSeq = 0.95
    Run(paf_dir, NbMatch, IdSeq)

    print('Terminée')