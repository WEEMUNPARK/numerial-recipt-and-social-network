#!/usr/bin/env python
# coding: utf-8

# In[6]:


import networkx as nx
import random
import pandas as pd
import numpy as np
data=pd.read_csv('marvel-unimodal-nodes.csv')
data
egdes=pd.read_csv('marvel-unimodal-edges.csv')





d={}
for i in range(len(data.Id)):
    d[data.Id[i]]=i+1

x=np.random.choice(list(d.values()),10)




p1=list(egdes.Source)#,egdes.Target))
p2=list(egdes.Target)
p3=list(egdes.Weight)
for i in range(len(p1)):
    p1[i]=d[p1[i]]
    p2[i]=d[p2[i]]
    p3[i]=p3[i]
new_edges=list(zip(p1,p2,p3))
new_edge=[]
for i in range(len(new_edges)):
    if new_edges[i][0] in x:
        if new_edges[i][1] in x:
            new_edge.append((new_edges[i][0],new_edges[i][1],0))
        else:
            new_edge.append((new_edges[i][0],new_edges[i][1],1))
    else:
        new_edge.append(new_edges[i])
        
        
        
new_edge=new_edge
G=nx.Graph()
for i in range(len(x)):
    G.add_node(x[i])
G.add_weighted_edges_from(new_edge)
pos = nx.spring_layout(G, k=0.1)
#plt.rcParams.update({'figure.figsize': (15, 10)})
nx.draw(G, alpha=0.05, with_labels=True)
#plt.show()



i=0
pas=[]
while i<100:
    a=np.random.choice(list(G.nodes()),2)
    
    pas.append(nx.shortest_path_length(G,a[0],a[1]))
    i+=1

np.mean(pas)


# In[ ]:




