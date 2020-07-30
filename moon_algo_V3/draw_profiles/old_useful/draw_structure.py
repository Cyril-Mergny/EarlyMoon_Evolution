#!/usr/bin/env python
# coding: utf-8

# In[2]:

import numpy as np  
import matplotlib.pyplot as plt
import pandas as pd


# In[12]:


# Import data
path = "../Outputs/"
filename = path + "profiles.txt"
names = ["int","r(m)","rho","vp","vs","eta"]
data = pd.read_csv(filename, skipinitialspace=True, sep=" ",index_col=False,names=names, skiprows=1)


# In[21]:


R = data["r(m)"]
R = [int(r/1e3) for r in R]

# In[34]:




fig, ax = plt.subplots(1,4,sharey=True, figsize=(12,5))
ax[0].plot(data["rho"], R)
ax[0].set_xlabel(r"Density $\rho$ (kg.m-3)")
ax[0].set_ylabel("Radius (km)")

ax[1].plot(data["vp"], R)
ax[1].set_xlabel("Vp (m.s-1)")
ax[2].plot(data["vs"], R)
ax[2].set_xlabel("Vs (m.s-1)")

ax[3].plot(data["eta"], R)
ax[3].set_xscale('log')
ax[3].set_xlabel("Viscosity $\eta$ (Pa.s)")
plt.savefig("../figures/structure_profiles.png")

#plt.show()



