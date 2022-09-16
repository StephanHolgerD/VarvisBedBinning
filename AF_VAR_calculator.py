import pandas as pd
from tqdm import tqdm
import numpy as np
from glob import glob
import yaml

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)




class AF_VAR_calculator():
    def __init__(self,YamlFile):
        self.config = YamlFile
        self.readConfig()    

    def readConfig(self):
        with open(self.config) as s:
            d=yaml.safe_load(s)
            self.BED =d['BED']
            self.chunksize= d['Targets']['chunksize']
            self.depth = d['Targets']['depth']
            self.af = d['Targets']['AF']
    
    
    
    def ChunkBed(self):
        chunk = self.chunksize
        bed = pd.read_csv(self.BED
                      ,sep='\t',names=['c','s','e'],usecols=[0,1,2])
        c_chr = list(bed['c'])
        c_start = list(bed['s'])
        c_end = list(bed['e'])
        c_chr = [c_chr[i:i+chunk] for i in range(0,len(c_chr),chunk)]
        c_end = [c_end[i:i+chunk] for i in range(0,len(c_end),chunk)]
        c_start = [c_start[i:i+chunk] for i in range(0,len(c_start),chunk)]

        self.c_chr,self.c_start,self.c_end=c_chr,c_start,c_end


  

    
    
    def BinVarAf(self,df):
        df = df[(df['reads']>=self.depth)&(df['alt_AF']>=self.af)]

        retL=[]
        test= np.array(df.values)
        startA = np.array(df['start'])
        endA = np.array(df['end'])
        chrA = np.array(df['chr'])
    
        for c,s,e in zip(self.c_chr,self.c_end,self.c_start):
            for cc in set(c):
                c_hits = [n for n,x in enumerate(c) if x==cc]
                ss=[s[i] for i in c_hits]
                ee=[e[i] for i in c_hits]
                cc=[c[i] for i in c_hits]
                min_s=ss[0]
                max_e=ee[-1]
                chrom=cc[0]
                i = np.where((startA>=min_s)&(endA<=max_e)&(chrA==chrom))
                i = test[i]
                afMean= np.average(i[:,3])            
                nVar = len(i)
                retL.append([chrom,min_s,max_e,nVar,afMean])
        return pd.DataFrame(retL,columns=['c','s','e','nVar','afM'])
    
    