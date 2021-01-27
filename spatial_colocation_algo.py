# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 17:02:58 2020

@author: HP
"""

from __future__ import division
import redis
import pandas as pd
from scipy import stats,integrate
import numpy as np
import re
import os
from  itertools import combinations
from operator import itemgetter
from itertools import combinations,groupby
from math import cos, asin, sqrt, pi

import time
start_time = time.time()

df=pd.read_csv("tree_data_PARK_AVE_SOUTH.csv")


df=df[["Botanical Name","Latitude","Longitude"]]
col_name='Botanical Name'
def distance(lat1, lon1, lat2, lon2):
    p = pi/180
    a = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p) * cos(lat2*p) * (1-cos((lon2-lon1)*p))/2
    return 12742 * asin(sqrt(a))

  
def remove(list_or_str): 
    if type(list_or_str)==list:
        pattern = '[0-9]'
        list_or_str = [re.sub(pattern, '', i) for i in list_or_str] 
        return list_or_str
    elif type(list_or_str)==str:
         list_or_str=''.join([i for i in list_or_str if not i.isdigit()])
         return list_or_str
  
def make_df_neighbors(crime,size):
    column_names=[]
    combination=list(combinations(crime, size))
    for i in range(len(combination)):
        combination[i]=list(combination[i])
        combination[i].sort()
        name=''
        for j in range(len(combination[i])):
            if (j==0):
                name=combination[i][j]
            else:
                name=name+"-"+combination[i][j]
        column_names.append(name)
    df_neighbors=pd.DataFrame(columns=column_names)
    return df_neighbors

def extract_event_types(df,col_name):
    crimes=list(df[col_name].unique())
    return crimes
    
def preprocess_df(df,col_name,crime):
    for i in crime:
        suffix=1
        for j in range(len(df)):
            if i==df[col_name][j]:
                df[col_name][j]=df[col_name][j]+str(suffix)
                suffix+=1
    df=df.rename(columns={col_name:"Event"})
    return df

def neighbors_size_2(df,neighbors_df, d_threshold):
    for i in range(len(df)):
        for j in range(i+1,len(df)):
            current_l=[]
            event_a=remove(df["Event"][i])
            event_b=remove(df["Event"][j])
            if(event_a!=event_b):
                current_l.append(df["Event"][i])
                if (distance(df["Latitude"][i],df["Longitude"][i],df["Latitude"][j],df["Longitude"][j])<=d_threshold):
                    current_l.append(df["Event"][j])
                    current_l.sort()
                    temp=[event_a,event_b]
                    temp.sort()
                    column_name=temp[0]+"-"+temp[1]
                    other_columns=list(set(neighbors_df.columns)-{column_name})
                    d=dict()
                    for k in list(neighbors_df.columns):
                        if k in other_columns:
                            d[k]=np.nan
                        else:
                            d[k]=current_l
                    neighbors_df = neighbors_df.append(d, ignore_index=True)
    neighbors_df = neighbors_df.apply(lambda x: pd.Series(x.dropna().values))
    return neighbors_df

def neighbors_size_k(size, neighbors_df, d_threshold, size_2_patterns, neighbors_df_size_2, 
                     neighbors_df_size_k_m_1, size_k_m_1_patterns):
    patterns_to_check=[]
    for i in range(len(size_k_m_1_patterns)):
        pattern1=[e for e in size_k_m_1_patterns[i]]
        pattern1_events=set(pattern1)
        for j in range(i+1,len(size_k_m_1_patterns)):
            pattern2=[e for e in size_k_m_1_patterns[j]]
            pattern2_events=set(pattern2)
            if len(pattern1_events.intersection(pattern2_events))!=0:
                temp=[]
                temp.append(size_k_m_1_patterns[i])
                temp.append(size_k_m_1_patterns[j])
                patterns_to_check.append(temp)
    if (len(patterns_to_check)!=0):
        for index in range(len(patterns_to_check)):
            size_k_m_1_patterns=patterns_to_check[index]
            for i in range(len(size_k_m_1_patterns)):
                pattern1=[e for e in size_k_m_1_patterns[i]]
                pattern1.sort()
                pattern1_name=("-").join(pattern1)
                pattern1_instances=list((neighbors_df_size_k_m_1[pattern1_name]).dropna())
                for j in range(i+1,len(size_k_m_1_patterns)):
                    pattern2=[e for e in size_k_m_1_patterns[j]]
                    pattern2.sort()
                    pattern2_name=("-").join(pattern2)
                    pattern2_instances=list((neighbors_df_size_k_m_1[pattern2_name]).dropna())
                    
                    for k in range(len(pattern1_instances)):
                        for l in range(len(pattern2_instances)):
                            x1=set(pattern1_instances[k])
                            x2=set(pattern2_instances[l])
                            common_elements=x1.intersection(x2)
                            if (len(common_elements)==size-2):
                                diff=list((x1-x2).union(x2-x1))
                                diff.sort()
                                diff_name=tuple(remove(diff))
                                diff_col_name=("-").join(remove(diff))
                                if diff_name in size_2_patterns:
                                    if (diff in [e for e in neighbors_df_size_2[diff_col_name]]):
                                        pattern=x1.union(x2)
                                        pattern_col_name=remove(("-").join(list(pattern)))
                                        other_columns=list(set(neighbors_df.columns)-{pattern_col_name})
                                        d=dict()
                                        for m in list(neighbors_df.columns):
                                            if m in other_columns:
                                                d[m]=np.nan
                                            else:
                                                d[m]=list(pattern)
                                        neighbors_df = neighbors_df.append(d, ignore_index=True)
                                else:
                                    continue
                            else:
                                continue
    neighbors_df = neighbors_df.apply(lambda x: pd.Series(x.dropna().values))
    return neighbors_df

 
def participation_index(df,neighbor):
    d_pi=dict()
    if (len(neighbor)!=0):
        for i in neighbor.columns:
            crimes=i.split("-")
            temp_l=list(neighbor[i].dropna())
            l_pr=[]
            for j in range(len(crimes)):
                joined_instances=[]
                for k in range(len(temp_l)):
                    temp=temp_l[k][j]
                    joined_instances.append(temp)
                joined_instances=set(joined_instances)
                crime_instances_joined=len(joined_instances)
                l=[]
            
                for m in range(len(df)):
                     if remove(df["Event"][m]) == crimes[j]:
                         l.append(df["Event"][m])
                for n in range(len(l)):
                    l[n] = int((int, re.findall(r'\d+', l[n]))[1][0])
                l.sort()
                total_crime_instances = l[-1]
                pr = crime_instances_joined/total_crime_instances
                l_pr.append(pr)
            d_pi[i] = min(l_pr)
    d_pi = dict(sorted(d_pi.items(), key = itemgetter(1),reverse = True))
    return d_pi

def spatial_colocation_miner_size2(df,crime, d_threshold,PARTICIPATION_THREOD):
    size=2
    neighbors_df=make_df_neighbors(crime,size) 
    neighbor=neighbors_size_2(df,neighbors_df, d_threshold)     
    p_i=participation_index(df,neighbor)
    l=list(p_i.items())
    output=[]
    for i in range(len(p_i)):
        current_output=[]
        if (l[i][1] >= PARTICIPATION_THREOD):
            current_output=(l[i][0]).split("-")
            output.append(tuple(current_output)) 
        else:
            break
    print("Size 2 done")
    return neighbor, output

def spatial_colocation_miner_sizek(df, crime, size, size_2_patterns, neighbors_df_size_2, 
                                   neighbors_df_size_k_m_1, size_k_m_1_patterns, 
                                   d_threshold, PARTICIPATION_THREOD):
    neighbors_df=make_df_neighbors(crime,size)
    neighbor=neighbors_size_k(size, neighbors_df, d_threshold, size_2_patterns, neighbors_df_size_2, 
                     neighbors_df_size_k_m_1, size_k_m_1_patterns)     
    p_i=participation_index(df,neighbor)
    l=list(p_i.items())
    output=[]
    for i in range(len(p_i)):
        current_output=[]
        if (l[i][1] >= PARTICIPATION_THREOD):
            current_output=(l[i][0]).split("-")
            output.append(tuple(current_output)) 
        else:
            break
    print("size", size, "done")    
    return neighbor, output



def spatial_colocation_miner(df, PARTICIPATION_THREOD = 0.4, d_threshold = 1):
    print(PARTICIPATION_THREOD, d_threshold)
    crime=extract_event_types(df,col_name)
    df=preprocess_df(df,col_name, crime)
    all_patterns=[] 
    global neighbors_df_size_2, size_2_patterns
    neighbors_df_size_2, size_2_patterns = spatial_colocation_miner_size2(df,crime, d_threshold, PARTICIPATION_THREOD)
    print(size_2_patterns)
    print("--- %s seconds ---" % (time.time() - start_time))
    neighbors_df_size_k_m_1, size_k_m_1_patterns = neighbors_df_size_2, size_2_patterns
    all_patterns.extend(size_k_m_1_patterns)
    size=3
    while (size_k_m_1_patterns!=[]):
        neighbors_df_size_k_m_1, size_k_m_1_patterns=spatial_colocation_miner_sizek(df, crime, size, 
                                                    size_2_patterns, neighbors_df_size_2, 
                                                    neighbors_df_size_k_m_1, size_k_m_1_patterns, d_threshold,
                                                    PARTICIPATION_THREOD)
        size+=1
        all_patterns.extend(size_k_m_1_patterns)
    return all_patterns
        
spatial_colocation_patterns=spatial_colocation_miner(df, PARTICIPATION_THREOD = 0.5, d_threshold = 0.015)
print("spatial_colocation_patterns: ",spatial_colocation_patterns)



print("--- %s seconds ---" % (time.time() - start_time))