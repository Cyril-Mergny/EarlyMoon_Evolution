#!/usr/bin/env python
# coding: utf-8

# In[41]:

#%matplotlib
import numpy as np  
import matplotlib.pyplot as plt
import pandas as pd


# In[42]:
path = "../Outputs/"

# Import data
filename = path + "Tcrust_profile.txt"
names = ["r(m)","T","dt"]
data = pd.read_csv(filename, skipinitialspace=True, sep=" ",index_col=False,names=names)

n = int(data["r(m)"][0]) #nbr of discrete radius in the crust
nbr_exp = int(data.shape[0]/n) # nbr of experiements

R = data["r(m)"][1:n+1] #Radius
R = [int(r/1e3) for r in R]


# In[43]:


#Extract Temperature and timestep
T_matrix = np.zeros((n,nbr_exp),dtype=float)
Time = []
for i in range(nbr_exp):
    T_matrix[:,i] = data["T"][1+i*n:(i+1)*n+1]
    Time.append(data["dt"][1+i*n])
Time = [time/365/86400/1e3 for time in Time]


# In[44]:
fig, ax = plt.subplots(figsize=(4,5))
for i in range(nbr_exp):
    ax.plot(T_matrix[:,i], R)
ax.set_ylabel("Radius (km)")
ax.set_xlabel(r"Temperature T(r) (K)")

plt.savefig("../figures/Tcrust_profile.png")
#plt.show()

