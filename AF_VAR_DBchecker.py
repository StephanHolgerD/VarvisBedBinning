import pandas as pd
from tqdm import tqdm
import numpy as np
from glob import glob
import yaml
import json

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)

class AF_VAR_DBchecker():
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
            self.DBfile = d['DB']


    def ReadDB(self):
        d= pd.read_json(self.DBfile)
        d=d.transpose()
        d['coords']=d.index
        d.index = list(range(len(d)))
        self.DB = d


    def checkPercentile(self,v,l):
        if None in l:
            return None
        if sum(l)==0:
            return  None
        for n,value in enumerate(l):
            if value>v:
                return n

        #perc = len([x for x in l if x<= v ])
        #return perc
    
    def Test_AF_VAR_df(self,df):
        df['key']=df['c'].astype(str) + '_' + df['s'].astype(str)+'_'+df['e'].astype(str)
        af_results=[]
        var_results=[]
        for k,af,var in zip(df['key'],df['afM'],df['nVar']):
            af=float(af)
            var = int(var)
            temp_d = self.DB[self.DB['coords']==k]
            
            af_distribution = temp_d['af'].values[0]
            var_distribution = temp_d['var'].values[0]
            percentile_af=self.checkPercentile(af,af_distribution)
            percentile_var=self.checkPercentile(var,var_distribution)
            af_results.append(percentile_af)
            var_results.append(percentile_var)
        df['Var_percentile']=var_results
        df['AF_percentile']=af_results
        return df