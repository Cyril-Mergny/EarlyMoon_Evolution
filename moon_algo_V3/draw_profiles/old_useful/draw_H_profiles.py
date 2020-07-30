#!/usr/bin/env python
# coding: utf-8

# In[401]:

import numpy as np  
import matplotlib.pyplot as plt
import pandas as pd


# Import data
filename = "../figures/H_profile.txt"
names = ["r(m)","Htide(W.m-3)","T(days)","n"]
data = pd.read_csv(filename, skipinitialspace=True, sep=" ", names=names, index_col=False)

n = data["n"][0] # nbr of radius
nbr_exp = int(data.shape[0]/n) # nbr of experiements

#Radius
R = data["r(m)"][:n]
R = [int(r/1e3) for r in R]

#Htide and Period
H_matrix = np.zeros((n,nbr_exp),dtype=float)
Periods = []
for i in range(nbr_exp):
    H_matrix[:,i] = data["Htide(W.m-3)"][i*n:(i+1)*n]
    Periods.append(round(data["T(days)"][i*n],2))

# # Individual Profile
#Plot individual profile
n_indiv = 0
T = Periods[n_indiv]
H = H_matrix[:,n_indiv]

# In[419]:
fig, (ax,ax2) = plt.subplots(1,2,figsize=(10,8))
ax.plot(H,R)
#ax.plot(H_matrix[:,0], R, label="Reference")
ax.set_xscale("log")
ax.set_title("H profile for T={} days".format(T))
ax.set_xlabel("Htide (W.m-3))", fontsize=12)
ax.set_ylabel("Radius (km)", fontsize=12)
#plt.savefig("../figures/Htide_profile.png")






# Import data
filename = "Hmu_mu.dat"
names = ["r(m)","Hmu_mu(W.m-3)","T(days)","n"]
data = pd.read_csv(filename, skipinitialspace=True, sep=" ", names=names, index_col=False)

n = data["n"][0]-1 # nbr of radius -1 here bc in K2_DT start at j=2
nbr_exp = int(data.shape[0]/n) # nbr of experiements

#Radius
R = data["r(m)"][:n]
R = [int(r/1e3) for r in R]

#Htide and Period
H_matrix = np.zeros((n,nbr_exp),dtype=float)
Periods = []
for i in range(nbr_exp):
    H_matrix[:,i] = data["Hmu_mu(W.m-3)"][i*n:(i+1)*n]
    Periods.append(round(data["T(days)"][i*n],2))


# # Individual Profile
#Plot individual profile
n_indiv = 0
T = Periods[n_indiv]
H = H_matrix[:,n_indiv]

# In[419]
ax2.plot(H,R)
#ax.plot(H_matrix[:,0], R, label="Reference")
ax2.set_xscale("log")
ax2.set_title("$H \mu \cdot \mu$ profile for T={} days".format(T))
#plt.legend()
ax2.set_xlabel("Htide (W.m-3))", fontsize=12)
ax2.set_ylabel("Radius (km)", fontsize=12)
plt.savefig("../figures/H_profiles.png")
#plt.show()






