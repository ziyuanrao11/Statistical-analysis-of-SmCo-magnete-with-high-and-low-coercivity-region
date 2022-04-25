# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:44:15 2019

@author: ThinkPad
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
# from comp_analysis_cluster import read_rrng, atom_filter
from tqdm import tqdm
import re
import numpy as np
# import matplotlib.pyplot as plt
import matplotlib.patches as patches
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import os
from mpl_toolkits import mplot3d
from scipy import stats
#import mpl_scatter_density
import seaborn as sns
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.stats import norm
import  pandas as pd    
#%%read the rrng file
def read_rrng(f):
    rf = open(f,'r').readlines()
    patterns = re.compile(r'Ion([0-9]+)=([A-Za-z0-9]+).*|Range([0-9]+)=(\d+.\d+) +(\d+.\d+) +Vol:(\d+.\d+) +([A-Za-z:0-9 ]+) +Color:([A-Z0-9]{6})')
    ions = []
    rrngs = []
    for line in rf:
        m = patterns.search(line)
        if m:
            if m.groups()[0] is not None:
                ions.append(m.groups()[:2])
            else:
                rrngs.append(m.groups()[2:])
    ions = pd.DataFrame(ions, columns=['number','name'])
    ions.set_index('number',inplace=True)
    rrngs = pd.DataFrame(rrngs, columns=['number','lower','upper','vol','comp','colour'])
    rrngs.set_index('number',inplace=True) 
    rrngs[['lower','upper','vol']] = rrngs[['lower','upper','vol']].astype(float)
    rrngs[['comp','colour']] = rrngs[['comp','colour']].astype(str)
    return ions, rrngs

def readpos(file_name):
    f = open(file_name, 'rb')
    dt_type = np.dtype({'names':['x', 'y', 'z', 'm'], 
                  'formats':['>f4', '>f4', '>f4', '>f4']})
    pos = np.fromfile(f, dt_type, -1)
    f.close()
    return pos

def label_ions(pos,rrngs):
    pos['comp'] = ''
    pos['colour'] = '#FFFFFF'
    pos['nature'] = ''
    count=0;
    for n,r in rrngs.iterrows():
        count= count+1;
        pos.loc[(pos.Da >= r.lower) & (pos.Da <= r.upper),['comp','colour', 'nature']] = [r['comp'],'#' + r['colour'],count]
    
    return pos



def atom_filter(x, Atom_range):
    Atom_total = pd.DataFrame()
    for i in range(len(Atom_range)):
        Atom = x[x['Da'].between(Atom_range['lower'][i], Atom_range['upper'][i], inclusive=True)]
        Atom_total = Atom_total.append(Atom)
        Count_Atom= len(Atom_total['Da'])   
    return Atom_total, Count_Atom  


#%%read the pos data to get the atom positions   
pos_file = 'B3 Low Hc R5076_53090. - complete.pos'
pos = readpos(pos_file)
dpos = pd.DataFrame({'x':pos['x'],
                            'y': pos['y'],
                            'z': pos['z'],
                            'Da': pos['m']})
#%%save the data   
dpos.to_csv('section_R5076_53090_LowHc.csv',index=False)  
del dpos
#%%cut the apt tip to 50 parts
# from tqdm import tqdm
all_data=pd.read_csv('section_R5076_53090_LowHc.csv')
sort_x = all_data.sort_values(by=['z'])    
max(sort_x['z'])
min(sort_x['z'])
sublength_x= int((max(sort_x['z'])-min(sort_x['z']))/50)
start = min(sort_x['z'])
end = min(sort_x['z'])+sublength_x
for i in tqdm(range(50)):
    # temp = sort_x.iloc[start:end]
    temp = sort_x[sort_x['z'].between(start, end, inclusive=True)]
    temp.to_csv('section_R5076_53090_LowHc/{}.csv'.format(i), index=False)
    start += sublength_x
    end += sublength_x
    print(end)    
    
    