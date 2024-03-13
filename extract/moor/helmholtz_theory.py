# Helmholtz bay theory plot
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

#import xarray as xr
import matplotlib.pyplot as plt
#import pandas as pd
import numpy as np

T_M2=12.42*3600 #M2 period
T_S2=12*3600 #S2 period

w_M2=(2*np.pi)/T_M2 #M2 frequency
w_S2=(2*np.pi)/T_S2 #S2 frequency
w=(w_M2+w_S2)/2 #forcing frequency

a=4e3*50 #cross-sectional area of inlet
A=8e3*40e3 #area of bay
g=9.81 #gravity
L_5=5e3
L_20=20e3
L_80=80e3

wo_5=np.sqrt((a*g)/(A*L_5))
wo_20=np.sqrt((a*g)/(A*L_20))
wo_80=np.sqrt((a*g)/(A*L_80))

wratio=np.linspace(0,2,200)

alpharatio_nofriction=1/(1-(wratio**2))

alpharatio_5=1/(1-((w/wo_5)**2))
alpharatio_20=1/(1-((w/wo_20)**2))
alpharatio_80=1/(1-((w/wo_80)**2))

Rworatio=0.4
alpharatio=1/np.sqrt(((1-(wratio**2))**2)+((Rworatio*wratio)**2))

plt.figure()
plt.plot(wratio,alpharatio_nofriction,label='R=0')
plt.plot(wratio,alpharatio,'--',label='R/wo=0.4')
plt.ylim(-3,3)
plt.xlim(0,2)
plt.plot([w/wo_5,w/wo_5],[-3,3],c='tab:red')
plt.plot([w/wo_20,w/wo_20],[-3,3],c='tab:green')
plt.plot([w/wo_80,w/wo_80],[-3,3],c='tab:purple')
plt.legend()

plt.text(w/wo_5,1.5,'5km',rotation='vertical')
plt.text(w/wo_20,1.5,'20km',rotation='vertical')
plt.text(w/wo_80,1.5,'80km',rotation='vertical')
plt.grid(True)
plt.show()




