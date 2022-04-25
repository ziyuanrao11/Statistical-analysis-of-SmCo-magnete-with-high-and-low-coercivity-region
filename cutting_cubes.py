# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:46:10 2019

@author: y.wei
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 17:31:55 2019

@author: y.wei
"""

import pandas as pd 
import numpy as np 
import os
from tqdm import tqdm
import re
import os
#folder = 'save'
#for filename in os.listdir(folder):
#    print(filename)
#    os.rename(folder+'/'+filename, folder+'/'+filename + '.csv')    
#s= pd.read_csv('save/chunk_0.csv')    
#y = x['comp'].dropna()
#y=y.unique()
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
    return ions,rrngs

def atom_filter(x, Atom_range):
    Atom_total = pd.DataFrame()
    for i in range(len(Atom_range)):
        Atom = x[x['Da'].between(Atom_range['lower'][i], Atom_range['upper'][i], inclusive=True)]
        Atom_total = Atom_total.append(Atom)
        Count_Atom= len(Atom_total['Da'])   
    return Atom_total, Count_Atom  
#Ola = sort_x.iloc[1:sublength_x]
#sort = Ola.sort_values(by=['z'])
#%%get the information of the elements
rrange_file = 'R5076_53078__SmCoFeCuZr_goodrangefile_20210309_Polin_WithoutHydrites_colors.rrng'
ions, rrngs = read_rrng(rrange_file)
Sm_range = rrngs[rrngs['comp']=='Sm:1']
Co_range = rrngs[rrngs['comp']=='Co:1']
Fe_range = rrngs[rrngs['comp']=='Fe:1']
Zr_range = rrngs[rrngs['comp']=='Zr:1']
Cu_range = rrngs[rrngs['comp']=='Cu:1']
#%%cut to cubes to be 2*2*2 nm
folder = 'section_R5076_53090_LowHc'
features = []
ratios = []
cube_size = 2
num=0
for filename in tqdm(os.listdir(folder)):
    print(filename)
    s= pd.read_csv(folder+'/'+filename)
    x_min = round(min(s['x']))
    x_max = round(max(s['x']))
    y_min = round(min(s['y']))
    y_max = round(max(s['y']))
    z_min = round(min(s['z']))
    z_max = round(max(s['z']))   
    p=[]
    x=[]
    folder_dir = 'section_R5076_53090_LowHc/cubic_{}'.format(num)  
    if not os.path.isdir(folder_dir):
        os.mkdir(folder_dir)
    num_cubic=0
    for i in range(z_min,z_max,cube_size):
        print(i)
        cubic = s[s['z'].between(i, i+cube_size, inclusive=True)]
        for j in range(y_min,y_max,cube_size):
            p = cubic[cubic['y'].between(j, j+cube_size, inclusive=True)]
    #        print(j)
            for k in range(x_min, x_max, cube_size):
#                print(k)
    #            print(k+5)
           
                x = p[p['x'].between(k,k+cube_size, inclusive=True)]
                if len(x['x'])>20:
                    name ='section_R5076_53090_LowHc/cubic_{}/{}.csv'.format(num,num_cubic) 
                    num_cubic+=1
                    x.to_csv(name, index=False)
    num+=1

                  
                    
#            print(x)
#    print('space')
#cubics_co = pd.DataFrame(cubics_co)
#cubics_co.columns= ['x','y','z']
#cubics_co.to_csv('cubics.csv',index=False)

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#plt.style.use('seaborn-white')
##plt.style.use('ggplot')  its nice!
#plt.rcParams['font.family'] = 'Arial'
#plt.rcParams['font.serif'] = 'Ubuntu'
#plt.rcParams['font.monospace'] = 'Ubuntu Mono'
#plt.rcParams['font.size'] = 18
#plt.rcParams['axes.labelsize'] = 18
#plt.rcParams['axes.titleweight'] = 'bold'
#plt.rcParams['axes.titlesize'] = 18
#plt.rcParams['xtick.labelsize'] = 18
#plt.rcParams['ytick.labelsize'] = 18
#plt.rcParams['legend.fontsize'] = 18
#plt.rcParams['figure.titlesize'] = 18





#fig = plt.figure(figsize=plt.figaspect(5)*2)
#ax = fig.gca(projection='3d')
##ax.axis('off')
##ax.axis('equal')
#ax.scatter(cubics_co['x'][:], cubics_co['y'][:], cubics_co['z'][:], c='r', s=5, alpha =0.8)
##ax.scatter(cubics_co['x_0'][:], cubics_co['y_0'][:], cubics_co['z_0'][:], c='r', s=5, alpha =0.8)
#
#ax.scatter(s['x'][::10000], s['y'][::10000], s['z'][::10000], c='b', s=5, alpha =0.4)

#
#Ti = zip(x[(x['comp']=='Ti:1')]['x'],x[(x['comp']=='Ti:1')]['y'],x[(x['comp']=='Ti:1')]['z'])
#Ti = list(Ti)
#Ti = [list(elem) for elem in Ti]
#points3 = np.array(Ti)
#
#mass_spect = s['Da']
#n, bins, patches = plt.hist(mass_spect, 10000, density=False, facecolor='r', alpha=0.5)
#plt.xlabel('X')
#plt.ylabel('Y')
##plt.tick_params(axis='x',which='minor',bottom=False)
#plt.title('Mass spectrum at grain boundary')
##plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(28, 38) 
##plt.minorticks_on()
#plt.yscale('log', nonposy='clip')
#plt.grid(linestyle='dotted')
#plt.show()

#cubics=[]
#for i in range(5):
#    print('i:')
#    print(i)
#    print('z')
#    cubic = total[total['z'].between(5*i, 5*i+5, inclusive=False)]
#    for j in range(3):
#        print('j:')
#        print(j)
#        cubic = cubic[cubic['y'].between(2*j, 2*j+2, inclusive=False)]
#        for k in range(3):
#            print('k:')
#            print(k)            
#            cubic = cubic[cubic['x'].between(2*k, 2*k+2, inclusive=False)]
#    print('space')


#plt.axis('off')   