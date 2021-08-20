#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:16:51 2021

FILTRADO MINUTAL

@author: inti
"""

#%%
import pandas as pd
import matplotlib.pyplot as plt
import pvlib as pv
import os

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

#%% INSPECCIONES VISUALES:

# #Serie temporal por separado
# plt.rcParams.update({'font.size': 12})
# fig, axs = plt.subplots(4, 1, figsize=(18,18))
# i=0

# for col in df1.columns[1:]:
#     axs[i].plot(df1.index,df0.loc[df0.index,col])
#     axs[i].set(ylabel=r'{}($W/m^2$)'.format(col),xlabel='Fecha',ylim=[-100,1750],title=col)
#     axs[i].grid()    
#     i = i+1
# plt.tight_layout()
# plt.savefig('{}/PX{:02d}_serieMinutal.png'.format(rutaFIG,PX))

# #Serie temporal GHI/DHI/GTI juntas
# df1[['GHI1','GHI2','DHI','GTI']].plot()
# plt.grid()
# plt.tight_layout()

# #fd vs kt pre procesado
# plt.rcParams.update({'font.size': 16})
# plt.figure(figsize=(10,8))
# plt.plot(kt,fd,'.',color='grey',label=r'$\theta_z > 80^{\circ}$')
# plt.plot(kt[df1.CZ>0.17],fd[df1.CZ>0.17],'.k',label=r'$\theta_z < 80^{\circ}$')
# plt.xlabel(r'$k_t$'); plt.ylabel(r'$f_d$')
# plt.xlim([-0.05, 1.2]); plt.ylim([-0.05, 1.2])
# plt.legend(markerscale=3,loc='center left', bbox_to_anchor=(1, 0.5))
# plt.grid()
# plt.tight_layout()
# plt.savefig('{}/PX{:02d}_fd_kt_prefilt.png'.format(rutaFIG,PX))


# #GHI y DHI vs CZ
plt.rcParams.update({'font.size': 16})
# fig, axs = plt.subplots(1, 3, figsize=(15,4))
# axs[0].plot(df1.CZ,df1.GHI1,'.k',markersize=2)
# axs[1].plot(df1.CZ,df1.GHI2,'.k',markersize=2)
# axs[2].plot(df1.CZ,df1.DHI,'.k',markersize=2)
# for i in range(0,3):
#     axs[i].set(ylabel=r'Irradiancia $(W/m^2)$',xlabel=r'$\cos \theta_z$',ylim=[-150,1750],title=df1.columns[i+1])
#     axs[i].grid()
# plt.tight_layout()
# plt.savefig('{}/PX{:02d}_GvsCZ.png'.format(rutaFIG,PX))


# #GHI1 vs GHI2
# plt.figure(figsize=(8,8))
# plt.plot(df1.GHI1,df1.GHI2,'.k',markersize=2)
# plt.plot([0,1500],[0,1500],'-r')
# plt.xlim([-100,1500]); plt.ylim([-100, 1500])
# plt.xlabel(r'$GHI_1(W/m^2)$');plt.ylabel(r'$GHI_2(W/m^2)$')
# plt.grid()
# plt.savefig('{}/PX{:02d}_GHI1vsGHI2.png'.format(rutaFIG,PX))
# plt.tight_layout()

#MASCARA DE DISCREPANCIA GHI1-GHI2
mskGHI12 = abs((df1.GHI1-df1.GHI2)/df1.GHI1) > 0.2

# plt.rcParams.update({'font.size': 16})
# plt.figure(figsize=(8,8))
# plt.plot(df1.GHI1,df1.GHI2,'.',color='grey',markersize=2,label='|GHI2-GHI1| < 0.1*GHI1')
# plt.plot(df1.GHI1[~mskGHI12],df1.GHI2[~mskGHI12],'.k',markersize=2,label='|GHI2-GHI1| > 0.1*GHI1')
# plt.plot([0,1500],[0,1500],'-r',label = 'GHI1 = GHI2')
# plt.xlim([-100,1500]); plt.ylim([-100, 1500])
# plt.xlabel(r'$GHI_1(W/m^2)$');plt.ylabel(r'$GHI_2(W/m^2)$')
# plt.legend(markerscale=4); plt.grid()
# plt.tight_layout()
# plt.savefig('{}/PX{:02d}_GHI1vsGHI2_msk.png'.format(rutaFIG,PX))

# plt.rcParams.update({'font.size': 16})
# plt.figure(figsize=(8,8))
# plt.plot(df1.GHI1,df1.GHI2,'.',color='grey',markersize=2,label=r'$|GHI_2-GHI_1| < 0.1 GHI_1$')
# plt.plot(df1.GHI1[~mskGHI12],df1.GHI2[~mskGHI12],'.k',markersize=2,label=r'$|GHI_2-GHI_1| > 0.1 GHI_1$')
# plt.plot([0,1500],[0,1500],'-r',label = 'GHI1 = GHI2')
# plt.xlim([-50,300]); plt.ylim([-50, 300])
# plt.xlabel(r'$GHI_1(W/m^2)$');plt.ylabel(r'$GHI_2(W/m^2)$')
# plt.legend(markerscale=4); plt.grid()
# plt.tight_layout()
# plt.savefig('{}/PX{:02d}_GHI1vsGHI2_msk_b.png'.format(rutaFIG,PX))

#fd vs kt
# plt.rcParams.update({'font.size': 16})
# plt.figure(figsize=(12,8))
# plt.plot(kt[df1.CZ>0.17],fd[df1.CZ>0.17],'.',color='grey',label=r'$|GHI_2-GHI_1| < 0.2 GHI_1$')
# plt.plot(kt[(df1.CZ>0.17)&(~mskGHI12)],fd[(df1.CZ>0.17)&(~mskGHI12)],'.k',label=r'$|GHI_2-GHI_1| > 0.1 GHI_1$')
# plt.xlabel(r'$k_t$'); plt.ylabel(r'$f_d$')
# plt.xlim([-0.05, 1.2]); plt.ylim([-0.05, 1.2])
# plt.legend(markerscale=3,loc='center left', bbox_to_anchor=(1, 0.5))
# plt.grid()
# plt.tight_layout()
# plt.savefig('{}/PX{:02d}_fd_kt_mskGHI12.png'.format(rutaFIG,PX))


#%% inspección de puntos malos en fd-kt

mskLOC = (fd > 0.15)&(fd < 0.23)&(kt > 0.25)&(kt < 0.40)
mskIdentifyPoints = (df1.CZ>0.17)&(~mskGHI12)&mskLOC

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(12,8))
plt.plot(kt[df1.CZ>0.17],fd[df1.CZ>0.17],'.',color='grey',label='datos')
plt.plot(kt[mskIdentifyPoints],fd[mskIdentifyPoints],'.r',label='dudosos')
plt.xlabel(r'$k_t$'); plt.ylabel(r'$f_d$')
plt.xlim([-0.05, 1.2]); plt.ylim([-0.05, 1.2])
plt.legend(markerscale=3,loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()
plt.tight_layout()

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8,8))
plt.plot(df1.GHI1,df1.GHI2,'.',color='grey',markersize=2,label=r'$|GHI_2-GHI_1| < 0.2 GHI_1$')
plt.plot(df1.GHI1[~mskGHI12],df1.GHI2[~mskGHI12],'.k',markersize=2,label=r'$|GHI_2-GHI_1| > 0.2 GHI_1$')
plt.plot(df1.GHI1[mskIdentifyPoints],df1.GHI2[mskIdentifyPoints],'.g',markersize=3,label='dudosos')
plt.plot([0,1500],[0,1500],'-r',label = 'GHI1 = GHI2')
plt.xlim([-100,1500]); plt.ylim([-100, 1500])
plt.xlabel(r'$GHI_1(W/m^2)$');plt.ylabel(r'$GHI_2(W/m^2)$')
plt.legend(markerscale=4); plt.grid()
plt.tight_layout()


#df1[['GHI1','GHI2','DHI']].plot()
plt.rcParams.update({'font.size': 12})
plt.figure()
plt.plot(df1.index,df1.GHI1,'-r',label='GHI1')
plt.plot(df1.index,df1.GHI2,'-b',label='GHI2')
plt.plot(df1.index,df1.DHI,'-',color='orange',label='DHI')
plt.plot(df1.loc[mskIdentifyPoints].index,df1[mskIdentifyPoints].GHI1,'*k', markersize=7,label='dudosos')
plt.legend(markerscale=3,loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()


plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(12,8))
plt.plot(kt[df1.CZ>0.17],fd[df1.CZ>0.17],'.',color='grey',label='datos')
plt.plot(kt.loc['2017-02-26'],fd.loc['2017-02-26'],'.r',label='26/02/2017')
plt.xlabel(r'$k_t$'); plt.ylabel(r'$f_d$')
plt.xlim([-0.05, 1.2]); plt.ylim([-0.05, 1.2])
plt.legend(markerscale=3,loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()
plt.tight_layout()
