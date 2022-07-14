# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 09:42:16 2022

@author: tnm12
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import sys
import Simulatorv2021 as RSDs #import simulator version 2.0.2.1

      
'''
##############################################################################
Simulations (standard rate constants)
##############################################################################
'''
ksd4 = 4e2/1e9

color = ['orange','aqua']

t_sim = np.linspace(0,6,1001)*3600 # seconds

k_txn = 0.0075 #trancription rate

con41 = [25,0]
IN4 = [50,0]
con51 = [0,25]
IN5 = [0,50]
REP_con = 500

model = RSDs.RSD_sim(5) # define the model instance and # of domains

# initialize species involved in the reaction
model.molecular_species('g{1,2}',DNA_con=25)
model.molecular_species('rep{2}',DNA_con=REP_con)

plt.subplot(2,4,1)
for n in range(len(IN4)):
    model.molecular_species('g{4,1}',DNA_con=con41[n])
    model.molecular_species('in{4}',DNA_con=IN4[n])
    model.molecular_species('g{5,1}',DNA_con=con51[n])
    model.molecular_species('in{5}',DNA_con=IN5[n])
    
    model.global_rate_constants(ktxn=k_txn) #globally changes transcription rates
    
    # simulate the model (input is simulation time)
    model.simulate(t_sim,smethod='LSODA') #noticed slight instability with RK45 in this example
    
    # pull out the species from the model solution to plot
    S2 = model.output_concentration('S{2}')
    
    fs = 12
    
    if n == 0:
        plt.plot(model.t/60,(S2/REP_con)*100,color=color[n],linewidth=2,linestyle='-')
    else:
        plt.plot(model.t/60,(S2/REP_con)*100,color=color[n],linewidth=2,linestyle='--')

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.ylim(-10,110)
plt.xlim(0,180)
ax1 = plt.gca()
ax1.xaxis.set_tick_params(which='both',size=3,width=1,direction='in',top='on')
#ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.yaxis.set_tick_params(which='both',size=3,width=1,direction='in',right='on')
plt.xlabel('Time (min)',fontsize=fs)
plt.ylabel('Reacted reporter (%)',fontsize=fs)
plt.title('standard rate constants')
plt.legend(['S2-G{4,1}v2','S2-G{5,1}v2'],frameon=False)

'''
##############################################################################
Simulations (slower I4 + 4_1 gate reaction)
##############################################################################
'''

t_sim = np.linspace(0,6,1001)*3600 # seconds

con41 = [25,0]
IN4 = [50,0]
con51 = [0,25]
IN5 = [0,50]

model = RSDs.RSD_sim(5) # define the model instance and # of domains



# specify DNA txn templates and reporters and concentrations

model.molecular_species('g{1,2}',DNA_con=25,krev=5/1e9)
model.molecular_species('R{2}',DNA_con=REP_con)

plt.subplot(2,4,3)
for n in range(len(IN4)):
    model.molecular_species('g{4,1}',DNA_con=con41[n],krsd=ksd4)
    model.molecular_species('in{4}',DNA_con=IN4[n])
    model.molecular_species('g{5,1}',DNA_con=con51[n])
    model.molecular_species('in{5}',DNA_con=IN5[n])
    
    model.global_rate_constants(ktxn=k_txn) #globally changes transcription rates
    
    # simulate the model (input is simulation time)
    model.simulate(t_sim,smethod='LSODA') #noticed slight instability with RK45 in this example
    
    # pull out the species from the model solution to plot
    S2 = model.output_concentration('S{2}')
    
    fs = 12
    
    plt.plot(model.t/60,(S2/REP_con)*100,color=color[n],linewidth=2,linestyle='--')

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.ylim(-10,110)
plt.xlim(0,180)
ax1 = plt.gca()
ax1.xaxis.set_tick_params(which='both',size=3,width=1,direction='in',top='on')
#ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.yaxis.set_tick_params(which='both',size=3,width=1,direction='in',right='on')
plt.xlabel('Time (min)',fontsize=fs)
plt.ylabel('Reacted reporter (%)',fontsize=fs)
plt.title('slower I4 + 4_1 gate reaction',fontsize=fs)
plt.legend(['S2-G{4,1}v2','S2-G{5,1}v2'],frameon=False)

plt.suptitle('SA22_SIfigure27B')