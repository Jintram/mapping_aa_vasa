#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 17:33:16 2020

@author: m.wehrens
"""

# Misc discarded code

#test1=pd.merge(df_list_sub[0],df_list_sub[1], how='inner')

#test1=df_list_sub[0].join(df_list_sub[1:], how='inner')




df_list_sub = np.transpose(np.squeeze(np.array([df.iloc[:,[3]].values for df in df_list])))
df = pd.DataFrame(data=df_list_sub, columns=colnames, index=)
B_test=df_list_sub[0]





list(df_list[0].index.values) 

all(df_list[0].iloc[:,0]==df_list[1].iloc[:,0])

df_list[1].iloc[:,3]
df_list[1].iloc[:,3].values