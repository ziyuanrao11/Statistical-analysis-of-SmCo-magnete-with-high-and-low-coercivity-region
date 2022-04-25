# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 18:58:54 2019

@author: y.wei
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 10:18:28 2019

@author: y.wei
"""

import pandas as pd
import re
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from numba import jit
from mpl_toolkits.mplot3d import Axes3D
from pandas import DataFrame
import seaborn as sns
import matplotlib
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
#%%define the elements in the data
rrange_file = 'R5076_53078__SmCoFeCuZr_goodrangefile_20210309_Polin_WithoutHydrites_colors.rrng'
ions, rrngs = read_rrng(rrange_file)
Sm_range = rrngs[rrngs['comp']=='Sm:1']
Co_range = rrngs[rrngs['comp']=='Co:1']
Fe_range = rrngs[rrngs['comp']=='Fe:1']
Zr_range = rrngs[rrngs['comp']=='Zr:1']
Cu_range = rrngs[rrngs['comp']=='Cu:1']
#%%read the 50 parts we cut
for i in tqdm(range(50)):
    folder = 'section_R5076_53090_LowHc/cubics_{}'.format(i)
    #N = []
    low_ratios = []
    high_ratios =[]
    low_total =pd.DataFrame()
    high_total = pd.DataFrame()
    for filename in os.listdir(folder):
    #    print(filename)
        x=[]
        x = pd.read_csv(folder+'/'+filename)
        Cu_total, Count_Cu = atom_filter(x, Cu_range)    
    #    Count_O = Count_ZrO+Count_HfO + Count_NbO + Count_TaO + Count_TiO
        N_x = len(x) 
    #    N_statistics = [Count_ZrO, Count_TiO, Count_Hf, Count_Nb, Count_Ta,
    #                    Count_Zr,Count_Ti, N_x]
    #    Ratio_statistics = [Count_ZrO/Count_total, Count_TiO/Count_total, Count_Hf/Count_total, Count_Nb/Count_total,
    #                        Count_Ta/Count_total, Count_Zr/Count_total, Count_Ti/Count_total]
        ratio = Count_Cu/N_x
    
    #    if ratio <= 0.01:
    #        low_D_ZrO_001 = low_D_ZrO_001.append(ZrO_total)  
    #        low_total=low_total.append(x)
    #        low_ratio = [Count_ZrO/N_x, Count_TiO/N_x, Count_Hf/N_x, Count_Nb/N_x,
    #                        Count_Ta/N_x, Count_Zr/N_x, Count_Ti/N_x]
    #        low_ratios.append(Ratio_statistics)
    #        print(ratio)
    #    if ratio >= 0.01:
    #        high_D_ZrO_001 = high_D_ZrO_001.append(x)
    #    if ratio >= 0.025:
    #        high_D_ZrO_003 = high_D_ZrO_003.append(ZrO_total)    
        if ratio >= 0.12:
    #        high_D_ZrO_004 = high_D_ZrO_004.append(O_total)
    #        temp= [ZrO_total, HfO_total, NbO_total, TaO_total, TiO_total]
    #        result = pd.concat(temp)
            high_total=high_total.append(x)
    high_total.to_csv('section_R5076_53090_LowHc/high_Cu_cubics_{}.csv'.format(i))
#%%save the data    
high_total = pd.DataFrame()    
#draw all of the plates
#%%plot 3d view
for i in tqdm(range(50)):   
    total=pd.read_csv('section_R5076_53090_LowHc/high_Cu_cubics_{}.csv'.format(i))
    high_total=high_total.append(total)
high_total = high_total.sort_values(by=['z'])  
high_total = high_total[::1000]
fig = plt.figure(figsize=[6.4,15])
ax = Axes3D(fig)
ax.scatter(high_total['x'], high_total['y'], high_total['z'], c='b', s=1, alpha =0.5)
plt.show()
#draw all of the data
all_data=pd.read_csv('SmCo.csv')
all_data = all_data.sort_values(by=['z'])  
all_data = all_data[::100]
fig = plt.figure(figsize=[6.4,15])
ax = Axes3D(fig)
ax.scatter(all_data['x'], all_data['y'], all_data['z'], c='b', s=1, alpha =0.5)
plt.show()
#count the total composition in the Cu-enriched phase
#%%calculate the composition in each voxle
Cu_Comp=[]
Sm_Comp=[]
Co_Comp=[]
Fe_Comp=[]
Zr_Comp=[]
for i in tqdm(range(50)):
    folder = 'section_R5076_53090_LowHc/cubics_{}'.format(i)
    for filename in os.listdir(folder):
        x=[]
        x = pd.read_csv(folder+'/'+filename)
        Cu_total, Count_Cu = atom_filter(x, Cu_range)   
        Sm_total, Count_Sm = atom_filter(x, Sm_range)  
        Co_total, Count_Co = atom_filter(x, Co_range)
        Fe_total, Count_Fe = atom_filter(x, Fe_range)  
        Zr_total, Count_Zr = atom_filter(x, Zr_range)
        
        N_x = len(x) 
        ratio_Cu = Count_Cu/N_x
        ratio_Sm = Count_Sm/N_x
        ratio_Co = Count_Co/N_x
        ratio_Fe = Count_Fe/N_x
        ratio_Zr = Count_Zr/N_x
        
        if ratio_Cu >= 0.12:
            Cu_Comp.append(ratio_Cu)
            Sm_Comp.append(ratio_Sm)
            Co_Comp.append(ratio_Co)
            Fe_Comp.append(ratio_Fe)
            Zr_Comp.append(ratio_Zr)
#%%save the data
Cu_Comp = DataFrame(Cu_Comp)
Cu_Comp.to_csv('section_R5076_53090_LowHc/1-5_phase_Cu_Comp.csv',index=False) 
Sm_Comp = DataFrame(Sm_Comp)
Sm_Comp.to_csv('section_R5076_53090_LowHc/1-5_phase_Sm_Comp.csv',index=False) 
Co_Comp = DataFrame(Co_Comp)
Co_Comp.to_csv('section_R5076_53090_LowHc/1-5_phase_Co_Comp.csv',index=False) 
Fe_Comp = DataFrame(Fe_Comp)
Fe_Comp.to_csv('section_R5076_53090_LowHc/1-5_phase_Fe_Comp.csv',index=False) 
Zr_Comp = DataFrame(Zr_Comp)
Zr_Comp.to_csv('section_R5076_53090_LowHc/1-5_phase_Zr_Comp.csv',index=False) 

info=[{'name':'Cu_Comp','color':'orange'},{'name':'Sm_Comp','color':'blue'},{'name':'Co_Comp','color':'red'},
      {'name':'Fe_Comp','color':'green'},{'name':'Zr_Comp','color':'purple'}]
TotalCompo = pd.DataFrame(columns=['Elements','Mean','Variance']) 
#%%plot the histgram
for i in range(len(info)): 
    Comp=pd.read_csv('section_R5076_53090_LowHc/1-5_phase_{}.csv'.format(info[i]['name']))
    Comp=Comp.to_numpy()
    Comp_mean = round(np.mean(Comp),3)
    Comp_std = round(np.std(Comp),3)
    print('the mean of {} is {}'.format(info[i]['name'],Comp_mean))
    print('the std of {} is {}'.format(info[i]['name'],Comp_std))
    TotalCompo=TotalCompo.append({'Elements':info[i]['name'],'Mean':Comp_mean,'Variance':Comp_std},ignore_index=True)
    plt.figure()
    sns.set_style('darkgrid')
    sns.distplot(Comp, bins=100, kde=True, hist_kws={'color':info[i]['color']},norm_hist=True)
    # weights = np.ones_like(Comp)/float(len(Comp))
    # plt.hist(Comp, bins=30, color =info[i]['color'],weights=weights)
    # sns.distplot(Cu_Comp, bins=30, kde=False, hist_kws={'color':'Purple'},norm_hist=True)
    # sns.distplot(Cu_Comp,hist=True,color='red',norm_hist=True,label="Cu_Comp")
    plt.xlabel ('{}'.format(info[i]['name']))
    plt.ylabel('probability')
    plt.title('Higstogram of {}'.format(info[i]['name']))
    # plt.xlim(0, 0.6)
    plt.grid(True)
    plt.show()
    plt.savefig('section_R5076_53090_LowHc/1-5_phase_{}.png'.format(info[i]['name']), format='png', dpi=300)
    
TotalCompo.to_csv('section_R5076_53090_LowHc/1-5_phase_total_Comp.csv',index=False)
plt.figure()
matplotlib.rcParams['font.sans-serif'] = ['SimHei']
matplotlib.rcParams['axes.unicode_minus'] = False
Compo = TotalCompo['Mean']
error_bar = TotalCompo['Variance']
error_params = dict(elinewidth=5, ecolor='steelblue', capsize=5)
plt.barh(range(5), Compo, height=0.7, color='steelblue', alpha=0.8,edgecolor='green', xerr=error_bar,error_kw=error_params, hatch='o')      # 从下往上画
plt.yticks(range(5), ['Cu', 'Sm', 'Co', 'Fe','Zr'])
plt.xlim(0,0.6)
plt.xlabel("Composition")
plt.title("Composition of different elements")
for x, y in enumerate(Compo):
    plt.text(y+0.01, x+0.1, '%s' % y)
plt.show()
plt.savefig('section_R5076_53090_LowHc/1-5_phase_total_composition.png', format='png', dpi=300)

#ax.scatter(ZrO_total['x'][::5], ZrO_total['y'][::5], ZrO_total['z'][::5], c='b', s=0.5, alpha =0.2)
#ax.scatter(low_D_ZrO_001['x'][::20], low_D_ZrO_001['y'][::20], 
#           low_D_ZrO_001['z'][::20],c='white', edgecolor='b', linewidth=0.5, alpha =0.5)

#ax.scatter(high_D_ZrO_003['x'], high_D_ZrO_003['y'], high_D_ZrO_003['z'],c='white', edgecolor='r', linewidth=0.5, alpha =0.5)
#ax.scatter(cubics['x'], cubics['y'], cubics['z'], c='k', s=0.5, alpha =0.2)
#high_total.to_csv('high_D_004_cubes.csv',index=False)

#from sklearn.cluster import MeanShift
#from scipy.spatial import ConvexHull
#
#
#Nb = zip(high_total['x'], high_total['y'], high_total['z'])
#Nb = list(Nb)
#Nb = [list(elem) for elem in Nb]
#points=np.array(Nb)
#mean_clusters = MeanShift(bandwidth=2).fit(points)
#mean_labels = mean_clusters.labels_
#n_clusters_mean = len(set(mean_labels)) - (1 if -1 in mean_labels else 0)
#mean_unique_labels = set(mean_labels)
#mean_colors = plt.cm.Spectral(np.linspace(0, 1, len(mean_unique_labels)))            
#mean_clusters =[]
#new_cluster = []
#for k, col in zip(mean_unique_labels, mean_colors):
#    if k != -1:
#        mean_clusters.append(points[mean_labels == k])  
#print('\n\n++ Meanshift Results')
#print('Estimated number of clusters:  Nb: %d' % (n_clusters_mean))
#
#for ii in range(0, 5, 5):
#    print('angle:{}'.format(ii))
#    fig = plt.figure(figsize=(10, 30))
#    hdb_axis1 = fig.add_subplot(111, projection='3d')
#    hdb_axis1.view_init(elev=100., azim=240)
#    hdb_axis1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#    hdb_axis1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#    hdb_axis1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#    hdb_axis1.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#    hdb_axis1.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#    hdb_axis1.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
##    hdb_axis1.axis('off')        
##    for i, col in zip(range(len(new_clusters)), new_colors):
##            print(col)
##            if np.mod(i,2) ==0:
##                hdb_axis1.scatter(new_clusters[i][:,0][::2], new_clusters[i][:,1][::2], 
##                                  new_clusters[i][:,2][::2], color=col, s= 1, alpha=0.5)
##        x = pd.DataFrame(new_clusters[i])
##        x.columns = ['X','Y','Z']
##        x.to_csv('clusters/new_cluster_{}.csv'.format(i), index=False)            
##    hdb_axis1.scatter(high_total['x'], high_total['y'], high_total['z'],
##                      color='b', linewidth=1, alpha =0.2)           
##
#    for i, col in zip(range(len(mean_clusters)), mean_colors):
#        if len(mean_clusters[i])>=4:
#            hull = ConvexHull(mean_clusters[i])
#            print(hull.volume)
#            if hull.volume>1:
#                new_cluster.append(mean_clusters[i])
#                for s in hull.simplices:
#                    hdb_axis1.plot(mean_clusters[i][s,0], mean_clusters[i][s,1], mean_clusters[i][s,2], linestyle='-', color= 'black', linewidth= 0.1)
#                hdb_axis1.scatter(mean_clusters[i][:,0], mean_clusters[i][:,1], 
#                                  mean_clusters[i][:,2], color=col, s=1, alpha=1)  #edgecolors='black', linewidths=2,  
#
#def is_p_inside_points_hull(cluster, p):
#    hull = ConvexHull(cluster)
#    points_inside=pd.DataFrame()
#    for i in tqdm(range(len(p))):
##            print(p.iloc[i])
#        p_co = p.iloc[i].drop(labels=['Da']) 
#        p_co = np.array(p_co).reshape(1,-1)
#        new_points = np.append(cluster, p_co, axis=0)
#        new_hull = ConvexHull(new_points)
#        if list(hull.vertices) == list(new_hull.vertices):
#           points_inside = points_inside.append(p.iloc[i])
##               print('point appended')
#    return points_inside
#
#
#total = pd.read_csv('total_lei.csv')
#for i, hdb_cluster in enumerate(new_cluster):
#     print(i)    
#     x = hdb_cluster
#     x_max = np.amax(x[:,0])
#     x_min = np.amin(x[:,0])
#     y_max = np.amax(x[:,1])
#     y_min = np.amin(x[:,1])
#     z_max = np.amax(x[:,2])
#     z_min = np.amin(x[:,2])
#     new_x = total[total['x'].between(x_min, x_max, inclusive=True)]
#     new_x = new_x[new_x['y'].between(y_min, y_max, inclusive=True)]
#     new_x = new_x[new_x['z'].between(z_min, z_max, inclusive=True)]
#     new_x1 = is_p_inside_points_hull(x, new_x)
#     name ='clusters/composition_cluster_{}.csv'.format(i)
#     print('cluster_composition_{}saved'.format(i))
#     new_x1.to_csv(name,index=False)
#import pickle
#
#with open('outfile', 'wb') as fp:
#    pickle.dump(new_cluster, fp)
#with open ('outfile', 'rb') as fp:
#    x = pickle.load(fp)    

    
    
#from hdbscan import HDBSCAN
#from scipy.spatial import ConvexHull
#Nb = zip(ZrO_total['x'], ZrO_total['y'], ZrO_total['z'])
#Nb = list(Nb)
#Nb = [list(elem) for elem in Nb]
#points=np.array(Nb)
#min_size=10
#n = 5
#hdb = HDBSCAN(min_cluster_size= min_size, cluster_selection_method='leaf',
#              allow_single_cluster=True, gen_min_span_tree=True).fit(points)
#hdb_labels = hdb.labels_
#n_clusters_hdb = len(set(hdb_labels)) - (1 if -1 in hdb_labels else 0)
#
##x = pd.read_csv('GB_datapoints.csv')
#
#print('\n\n++ HDBSCAN Results')
#print('Estimated number of clusters:  Nb: %d' % (n_clusters_hdb))
#hdb_unique_labels = set(hdb_labels)
#hdb_colors = plt.cm.hsv(np.linspace(0, 1, len(hdb_unique_labels)))
#hdb_clusters =[]
#for k, col in zip(hdb_unique_labels, hdb_colors):
#    if k != -1:
#        hdb_clusters.append(points[hdb_labels == k])      
#for ii in range(0, 5, 5):
#    print('angle:{}'.format(ii))
#    fig = plt.figure(figsize=(10, 30))
#    hdb_axis1 = fig.add_subplot(111, projection='3d')
#    hdb_axis1.view_init(elev=0., azim=240)
#    hdb_axis1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#    hdb_axis1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#    hdb_axis1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#    hdb_axis1.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#    hdb_axis1.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#    hdb_axis1.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
##    hdb_axis1.axis('off')        
##    for i, col in zip(range(len(new_clusters)), new_colors):
##            print(col)
##            if np.mod(i,2) ==0:
##                hdb_axis1.scatter(new_clusters[i][:,0][::2], new_clusters[i][:,1][::2], 
##                                  new_clusters[i][:,2][::2], color=col, s= 1, alpha=0.5)
##        x = pd.DataFrame(new_clusters[i])
##        x.columns = ['X','Y','Z']
##        x.to_csv('clusters/new_cluster_{}.csv'.format(i), index=False)            
#    hdb_axis1.scatter(points[:,0], points[:,1], points[:,2], color='black', s= 3, alpha=0.5)            
#
#    for i, col in zip(range(len(hdb_clusters)), hdb_colors):
#        hull = ConvexHull(hdb_clusters[i])
##        for s in hull.simplices:
##            hdb_axis1.plot(hdb_clusters[i][s,0], hdb_clusters[i][s,1], hdb_clusters[i][s,2],
##                           linestyle='-', color= 'black', linewidth= 0.2)    
#        hdb_axis1.scatter(hdb_clusters[i][:,0], hdb_clusters[i][:,1], 
#                          hdb_clusters[i][:,2], color='white', edgecolors=col, linewidths=1, s=5, alpha=0.5)  #edgecolors='black', linewidths=2,      
##        