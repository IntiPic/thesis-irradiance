#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:16:51 2021

FILTRADO MINUTAL

@author: inti
"""

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pvlib as pv
import datetime as dt
import os
from fun import filtros as ms

# %% 1 - CARGA Y PREPARACIÓN DE DATOS
est = 'AR'
PX = 1
LATdeg = -30.39820 ; LONdeg = -56.51210; Zsnm = 136 #AR
tz = 'Etc/GMT+3'

loc = pv.location.Location(LATdeg,LONdeg,tz,Zsnm,est)

#%% Acceso a datos
cwd = os.getcwd()
rutaPX = cwd.replace('scripts/{}'.format(est), 'datos/{}/PX{:02d}'.format(est,PX)) 
rutaFIG = cwd.replace('scripts/{}'.format(est), 'fig/{}'.format(est))
rutaOUT = cwd.replace('scripts/{}'.format(est), 'out/{}'.format(est)) 

fileName = '{}/PX{:02d}_{}_2017_2018.csv'.format(rutaPX,PX,est)

#%% #### Levanto datos en un dataFrame
df0 = pd.read_csv(fileName, index_col = 'Fecha')
df0.index = pd.to_datetime(df0.index)

# Copia y ceros en la noche
df1 = df0.copy(deep=True)
df1[df0.CZ<0] = 0

#Aplico la calibración
kg = 0.997/1.006 
kd = 1.079/1.048
df1['GHI2'] = kg*df0['GHI2'] + 10.603
df1['DHI'] = kd*df0['DHI'] + 10.212

#%% #### Calculo fd y kt, y antes la DNI extraterrestre 
Iob = pv.irradiance.get_extra_radiation(df0.index.dayofyear, solar_constant=1367)
Iob = pd.Series(Iob,index=df1.index)

fd = df1.DHI / df1.GHI2
kt = df1.GHI1 / (Iob * df1.CZ)

#%% Determinación parámetros BSRN

#%% ### MASCARAS NAN
mskGHInan = (~df1.GHI1.isna()) & (~df1.GHI2.isna()) 
mskDHInan = ~df1.DHI.isna()

mskNAN_PH = (mskGHInan & mskDHInan)
# BASE = DIURNOS + No NAN
mskBASE_PH = (df1.CZ>0)&mskNAN_PH

a1 = 1.2; b1=1.2; c1=30
a2 = 0.75; b2=1.2; c2=30

CZaux = np.arange(0,1.05,0.05)

#%% PROCEDIMIENTO PARA HALLAR PARÁMETROS BSRN
'''
#-------------------------------------------------------------------------------------------#
                  GHI 
Fijo b = 1, varío a.  Grafico contra datos, elijo a
#-------------------------------------------------------------------------------------------#
'''

A1 = np.arange(0.9,1.5,0.1); B = np.arange(1.0,1.5,0.1)
b = 1.0
plt.figure(); plt.title('GHI = Csc*a*CZ^b + {} / a variable, b = 1'.format(c1))
plt.plot(df1.CZ[mskBASE_PH],df1.GHI2[mskBASE_PH],'.',color='grey', label = 'datos', markersize=2)
plt.ylim([-100,1750])
plt.xlabel('CZ'); plt.ylabel('Irradiancia (W/m^2)')
for a in A1:
    GHImax_alt = 1367*a*(CZaux**b) + c1
    plt.plot(CZaux,GHImax_alt,'-',label='a = {:.1f}'.format(a), markersize=2.0)

plt.legend(markerscale=4,loc='best'); plt.grid()

# Fijo a, varío b. Grafico contra datos, elijo b
a = 1.10
plt.figure(); plt.title('GHI = Csc*{:.2f}*CZ^b + {} / b variable'.format(a,c1))
plt.plot(df1.CZ[mskBASE_PH],df1.GHI2[mskBASE_PH],'.',color='grey', label = 'datos', markersize=2)
plt.ylim([-100,1750])
plt.xlabel('CZ'); plt.ylabel('Irradiancia (W/m^2)')
for b in B:
    GHImax_alt = 1367*a*(CZaux**b) + c1
    plt.plot(CZaux,GHImax_alt,'-',label='b = {:.1f}'.format(b), markersize=2.0)

plt.legend(markerscale=4,loc='best'); plt.grid()

#-------------------------------------------------------------------------------------------#
#                  DHI
# Fijo b = 1, varío a.  Grafico contra datos, elijo a
A2 = np.arange(0.6,1.3,0.1); B = np.arange(1.0,1.5,0.1)
b=1.0
plt.figure(); plt.title('DHI = Csc*a*CZ^b + {}/ a variable, b = 1'.format(c2))
plt.plot(df1.CZ[mskBASE_PH],df1.DHI[mskBASE_PH],'.',color='grey', label = 'datos', markersize=2)
plt.ylim([-100,1750])
plt.xlabel('CZ'); plt.ylabel('Irradiancia (W/m^2)')
for a in A2:
    GHImax_alt = 1367*a*(CZaux**b) + c2
    plt.plot(CZaux,GHImax_alt,'-',label='a = {:.1f}'.format(a), markersize=2.0)

plt.legend(markerscale=4,loc='best'); plt.grid()

# # Fijo a, varío b. Grafico contra datos, elijo b
a = 0.6
plt.figure(); plt.title('DHI = Csc*{:.2f}*CZ^b + {} / b variable'.format(a,c2))
plt.plot(df1.CZ[mskBASE_PH],df1.DHI[mskBASE_PH],'.',color='grey', label = 'datos', markersize=2)
plt.ylim([-100,1000])
plt.xlabel('CZ'); plt.ylabel('Irradiancia (W/m^2)')
for b in B:
    GHImax_alt = 1367*a*(CZaux**b) + c2
    plt.plot(CZaux,GHImax_alt,'-',label='b = {:.1f}'.format(b), markersize=2.0)

plt.legend(markerscale=4,loc='best'); plt.grid()


#%% Filtrado BSRN GHI-DHI final

mskF0 = ms.F0(df1.GHI1)
mskF02 = ms.F0(df1.GHI2)
mskF1 = ms.F1(df1.DHI)

a1 = 1.15; b1=1.2; c1=30
a2 = 0.6; b2=1.1; c2=30
GHImax = Iob*a1*(df1.CZ**b1) + c1
DHImax = Iob*a2*(df1.CZ**b2) + c2
CZaux = np.arange(0,1.05,0.05)
GHImax_alt = 1367*a1*(CZaux**b1) + c1
DHImax_alt = 1367*a2*(CZaux**b2) + c2

mskF3 = ms.F3(df1.GHI1,df1.CZ,Iob,a1,b1,c1) #EXTREMELY RARE
mskF32 = ms.F3(df1.GHI2,df1.CZ,Iob,a1,b1,c1) #EXTREMELY RARE
mskF4 = ms.F4(df1.DHI,df1.CZ,Iob,a2,b2,c2) #EXTREMELY RARE

mskGHI1 = mskF0|mskF3
mskGHI2 = mskF02|mskF32
mskDHI = mskF1|mskF4

plt.rcParams.update({'font.size': 16})
#GHI y DHI vs CZ
fig, axs = plt.subplots(1, 3, figsize=(18,5))
axs[0].plot(df1.CZ,df1.GHI1,'.',color='grey',markersize=2,label='pasan BSRN')
axs[0].plot(df1[mskGHI1].CZ,df1[mskGHI1].GHI1,'.r',markersize=2, label= 'descarte BSRN')
axs[0].set(ylabel=r'Irradiancia $(W/m^2)$',xlabel=r'$\cos\theta_z)$',ylim=[-100,1750],title='GHI1')
axs[0].grid(); axs[0].legend(markerscale=4)
axs[1].plot(df1.CZ,df1.GHI2,'.',color='grey',markersize=2,label='pasan BSRN')
axs[1].plot(df1[mskGHI2].CZ,df1[mskGHI2].GHI2,'.r',markersize=2, label= 'descarte BSRN')
axs[1].set(ylabel=r'Irradiancia $(W/m^2)$',xlabel=r'$\cos\theta_z)$',ylim=[-100,1750],title='GHI2')
axs[1].grid(); axs[1].legend(markerscale=4)
axs[2].plot(df1.CZ,df1.DHI,'.',color='grey',markersize=2,label='pasan BSRN')
axs[2].plot(df1[mskDHI].CZ,df1[mskDHI].DHI,'.r',markersize=2,label='descarte BSRN')
axs[2].set(ylabel=r'Irradiancia $(W/m^2)$',xlabel=r'$\cos\theta_z)$',ylim=[-100,1750],title='DHI')
axs[2].grid(); axs[2].legend(markerscale=4)
plt.tight_layout()
#plt.savefig('{}/PX{:02d}_filtradoBSRN.png'.format(rutaFIG,PX))