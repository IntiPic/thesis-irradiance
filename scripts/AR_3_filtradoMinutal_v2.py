#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import pvlib as pv
import matplotlib.pyplot as plt
from fun import modelos as md
from fun import filtros as ms

#%% ### 1 - CARGA Y PREPARACIÓN DE DATOS
est = 'AR'
PX = 1
LATdeg = -30.39820 ; LONdeg = -56.51210; Zsnm = 136 #AR
BETAdeg = 30
tz = 'Etc/GMT+3'

loc = pv.location.Location(LATdeg,LONdeg,tz,Zsnm,est)
#%% #### Acceso a datos
cwd = os.getcwd()
rutaPX = cwd.replace('scripts/{}'.format(est),'') + 'datos/{}/PX{:02d}/'.format(est,PX) 
rutaMSK = cwd.replace('scripts/{}'.format(est), 'msk/{}'.format(est))
rutaFIG = cwd.replace('scripts/{}'.format(est), 'fig/{}'.format(est))
rutaOUT = cwd.replace('scripts/{}'.format(est), 'out/{}'.format(est)) 

fileName = rutaPX + 'PX{:02d}_{}_2017_2018.csv'.format(PX,est)
#%% #### Levanto datos en un dataFrame
df0 = pd.read_csv(fileName, index_col = 'Fecha')
df0.index = pd.to_datetime(df0.index)

#%% Copia y ceros en la noche
df1 = df0.copy(deep=True)
df1[df0.CZ<0] = 0

#Aplico la calibración
kg = 0.997/1.006 
kd = 1.079/1.048
df1['GHI2'] = kg*df0['GHI2'] + 10.603
df1['DHI'] = kd*df0['DHI'] + 10.212



#%% #### Calculo fd y kt, y antes la DNI extraterrestre 
Iob = pv.irradiance.get_extra_radiation(df1.index.dayofyear, solar_constant=1367)
Iob = pd.Series(Iob,index=df1.index)

fd = df1.DHI / df1.GHI2
kt = df1.GHI1 / (Iob * df1.CZ)

#%% Calculo las componentes de cielo claro y su GTIcsk para posterior filtrado
TL_01 = 1.5; TL_02 = 2.5
ciclosTL = pd.read_csv(rutaPX + '/TL_ciclos_por_estacion_esra', sep=",",index_col='MES')

TLcic = md.getTL(ciclosTL,est,df1.index.dayofyear)
TLcic = pd.Series(data=TLcic, index=df1.index)


alt_solar = np.rad2deg(np.arcsin(df1.CZ))

GHIcs1,_,DNIcs1 = md.generaComponentesESRA(alt_solar,Iob,Zsnm,TL=TLcic-TL_01) #Más estricto, para F6 y F16
GHIcs2,DHIcs2,DNIcs2 = md.generaComponentesESRA(alt_solar,Iob,Zsnm,TL=TLcic-TL_02) #Mas laxo, para F19

#GHIsum = df1.DNI1*df1.CZ + df1.DHI

AM = md.am_kastenYoung(Zsnm,alt_solar)

#%%----------------------------------------------------------------------------
### MASCARAS
#---------------------------------------------------------------------------
# fechaCambioInclinacion1 = '2018/03/13 23:50:00'
# fechaCambioInclinacion2 = '2018/03/14 00:00:00'
# mskCambioInc1 = (df1.index <= fechaCambioInclinacion1)
# mskCambioInc2 = (df1.index >= fechaCambioInclinacion2)


#%%Calculo la GTIcs haciendo ESRA + ISO 
solar_pos_1s = loc.get_solarposition(df1.index)

theta = pv.irradiance.aoi(BETAdeg,0,solar_pos_1s['zenith'], solar_pos_1s['azimuth'])


GTIcs2 = (np.maximum(DNIcs2 * np.cos(np.deg2rad(theta)),0)
            + DHIcs2*0.5*(1+np.cos(np.deg2rad(BETAdeg))) 
                + GHIcs2*0.25*0.5*(1-np.cos(np.deg2rad(BETAdeg))))



#%% ### MASCARAS NAN
mskGHInan = (~df1.GHI1.isna()) & (~df1.GHI2.isna()) 
mskDHInan = ~df1.DHI.isna()
mskGTInan = ~df1.GTI.isna()
mskNAN_PH = (mskGHInan & mskDHInan)
mskNAN_PI = mskNAN_PH & mskGTInan

#%%Genero una máscara de fechas para probar la distribución estacional
t1_30 = '2016/01/01 00:00:00'; t2_30 =  '2018/02/01 23:50:00'
mskFECHAS_30 = (t1_30 <= df1.index) & (df1.index <= t2_30) 


#%%
# ### APLICO FILTROS

# BASE = DIURNOS + No NAN
mskBASE_PH = (df1.CZ>0)&mskNAN_PH
mskBASE_PI = (df1.CZ>0)&mskNAN_PI

mskF0 = ms.F0(df1.GHI1,0); mskF0_2 = ms.F0(df1.GHI2,0)
mskF1 = ms.F1(df1.DHI,0)
# mskF2 = ms.F2(df1.DNI1)

mskF3 = ms.F3(df1.GHI1,df1.CZ,Iob,1.15,1.2,30) #EXTREMELY RARE
mskF3_2 = ms.F3(df1.GHI2,df1.CZ,Iob,1.15,1.2,30) #EXTREMELY RARE

mskF4 = ms.F4(df1.DHI,df1.CZ,Iob,0.60,1.1,30) #EXTREMELY RARE
    
# mskF5 = ms.F5(df1.DNI1,Iob)
# mskF6 = ms.F6(df1.DNI1,DNIcs1)

# #clausura
# mskF7 = ms.F7(df1.GHI1,GHIsum,df1.CZ,0.08,-50) #último argumento es el Gsum minimo
# mskF8 = ms.F8(df1.GHI1,GHIsum,df1.CZ,0.15,-50) #último argumento es el Gsum minimo

#fd min y max
mskF9 = ms.F9(df1.GHI2,df1.DHI,df1.CZ,1.05,-10) #último argumento es el GHI minimo
mskF10 = ms.F10(df1.GHI2,df1.DHI,df1.CZ,1.10,-10) #último argumento es el GHI minimo

# mskF11 = ms.F11(df1.GHI1,df1.DHI,df1.DNI1,df1.CZ,50)

mskF12 = ms.F12(df1.DHI,850)
#mskF13 = ms.F13(df1.DHI,Iob,df1.CZ,0.6) # Redundante con F3 (otros parámetros)

#Rectángulo
mskF14 = ms.F14(df1.GHI2,df1.DHI,df1.CZ,Iob,0.15,0.90) #max (kt,fd) - Elimina taco
mskF15 = ms.F15(df1.GHI2,df1.DHI,df1.CZ,Iob,0.60,0.98) #min (kt,fd) - Elimina fd=1 kt alto

# mskF16 = ms.F16(df1.DHI,GHIsum,GHIcs1,-50,0.80) #Desalineacion (keep an eye)

mskF17 = ms.F17(df1.GHI1,Iob,df1.CZ,AM,1.05) # Usa Kt perez
mskF17_2 = ms.F17(df1.GHI2,Iob,df1.CZ,AM,1.05) # Usa Kt perez

mskF18 = ms.F18(df1.GHI1,df1.DHI,DNIcs1,df1.CZ)
mskF18_2 = ms.F18(df1.GHI2,df1.DHI,DNIcs1,df1.CZ)

#Filtros GTI 
mskF19 = ((df1.GTI > GTIcs2) & (df1.GTI > 600))  # Filtro de GTI
mskF20 = df1.GTI < -4  

#Filtros sol bajo PH y PI
mskF21 = np.rad2deg(np.arccos(df1.CZ)) > 80
mskF22 = theta > 80

mskGHI12 = abs((df1.GHI1-df1.GHI2)/df1.GHI1) > 0.2

# Filtro de fechas donde pasan cosas raras (eclipse, etc)

# f11 = '2020-12-14 00:00:00'; f12 = '2020-12-14 23:59:00' #Eclipse
# f21 = '2015-10-29 00:00:00'; f22 = '2015-10-29 23:59:00' #fd = 1 y kt alto no filtrado por F16
# f31 = '2015-06-11 00:00:00'; f32 = '2015-06-14 23:59:00' #pico en fd=0.1 y kt = 0.8

# mskF23 = (((df1.index >= f11) & (df1.index <= f12)) | 
#           ((df1.index >= f21) & (df1.index <= f22)) |
#           ((df1.index >= f31) & (df1.index <= f32)))

#%% Creación de df de máscaras

mskDict = ({'mskBASE_PH':mskBASE_PH,'F0':mskF21,'F1':mskF0|mskF0_2,'F2':mskF1,
            'F4':mskF3|mskF3_2,'F5':mskF4,
            'F10':mskF9,'F11':mskF10,'F13':mskF12,'F14':mskF14,
            'F15':mskF15,'F17':mskF17|mskF17_2,'F18':mskF18|mskF18_2,'F19':mskF19,
            'F20':mskF20,'F21':mskF22,'Fe':mskGHI12})


## #### Guardo los resultados del filtrado en un dataframe y lo salvo en la carpeta "msk"
#

mascarasResult = pd.DataFrame(index = df1.index, columns = mskDict.keys())

for msk in mskDict.keys():
    mascarasResult[msk] = mskDict[msk]

#%% GUARDO SALIDAS
try:
    os.mkdir(rutaMSK)
except: pass
# try: 
#     os.mkdir('{}/{}/'.format(rutaMSK,est))
# except: pass

mascarasResult01 = mascarasResult.copy(deep=True)

mascarasResult = mascarasResult.astype(int)
mascarasResult.to_csv('{}/PX{:02d}_{}_MSK.csv'.format(rutaMSK,PX,est))

#%% Primeras pruebas gráficas
# plt.rcParams.update({'font.size': 12})
# fig, axs = plt.subplots(1, 3, figsize=(18,4))
# axs[0].plot(df1.CZ[mskBASE_PH],df1.GHI1[mskBASE_PH],'.',color='grey', label='descarte', markersize=2)
# axs[0].plot(df1.CZ[mskBASE_PH&(~mskF0)&(~mskF3)],df1.GHI1[mskBASE_PH&(~mskF0)&(~mskF3)],'.k',label='pasan', markersize=2)
# axs[0].plot(df1.CZ,1367*1.2*df1.CZ**1.1 + 40,'.r',label='Gcs*a*CZ^b+c', markersize=2)
# axs[0].set(xlabel='GHI(W/m^2)',ylabel='cos(theta_z)',title='F0 y F3',ylim=[-100,1750])
# axs[0].legend(markerscale=4); axs[0].grid()
# axs[1].plot(df1.CZ[mskBASE_PH],df1.DHI[mskBASE_PH],'.',color='grey', label='descarte', markersize=2)
# axs[1].plot(df1.CZ[mskBASE_PH&(~mskF1)&(~mskF4)],df1.DHI[mskBASE_PH&(~mskF1)&(~mskF4)],'.k',label='pasan', markersize=2)
# axs[1].plot(df1.CZ,1367*0.75*df1.CZ**1.2 + 30,'.r',label='Gcs*a*CZ^b+c', markersize=2)
# axs[1].set(xlabel='DHI(W/m^2)',ylabel='cos(theta_z)',title='F1 y F4',ylim=[-100,1750])
# axs[1].legend(markerscale=4); axs[1].grid()
# axs[2].plot(df1.CZ[mskBASE_PH],df1.DNI1[mskBASE_PH],'.',color='grey', label='descarte', markersize=2)
# axs[2].plot(df1.CZ[mskBASE_PH&(~mskF2)&(~mskF5)&(~mskF6)],df1.DNI1[mskBASE_PH&(~mskF2)&(~mskF5)&(~mskF6)],'.k',label='pasan', markersize=2)
# axs[2].set(xlabel='DNI(W/m^2)',ylabel='cos(theta_z)',title='F2, F5 y F6',ylim=[-250,1750])
# axs[2].legend(markerscale=4); axs[2].grid()
# plt.tight_layout()

#%% FILTROS por variable 

flagGHI = ['F0','F1','F4','F10','F11','F14','F15','F17','F18','Fe']
flagDHI = ['F0','F2','F5','F10','F11','F13','F14','F15','F18','Fe']


mskGHI = (mascarasResult[flagGHI]).any(axis='columns')
mskDHI = (mascarasResult[flagDHI]).any(axis='columns')

# plt.rcParams.update({'font.size': 14})
# fig, axs = plt.subplots(1, 3, figsize=(18,4))
# axs[0].plot(df1.CZ[mskBASE_PH],df1.GHI1[mskBASE_PH],'.',color='red', label='descarte', markersize=2)
# axs[0].plot(df1.CZ[mskBASE_PH&(~mskGHI)],df1.GHI1[mskBASE_PH&(~mskGHI)],'.b',label='pasan', markersize=2)
# axs[0].set(ylabel=r'$GHI_1(W/m^2)$',xlabel=r'$\cos(\theta_z)$',title=r'$GHI_1$',ylim=[-150,1500])
# axs[0].legend(markerscale=4, loc='best'); axs[0].grid()
# axs[1].plot(df1.CZ[mskBASE_PH],df1.GHI2[mskBASE_PH],'.',color='red', label='descarte', markersize=2)
# axs[1].plot(df1.CZ[mskBASE_PH&(~mskGHI)],df1.GHI2[mskBASE_PH&(~mskGHI)],'.b',label='pasan', markersize=2)
# axs[0].set(ylabel=r'$GHI_2(W/m^2)$',xlabel=r'$\cos(\theta_z)$',title=r'$GHI_2$',ylim=[-150,1500])
# axs[1].legend(markerscale=4, loc='best'); axs[1].grid()
# axs[2].plot(df1.CZ[mskBASE_PH],df1.DHI[mskBASE_PH],'.',color='red', label='descarte', markersize=2)
# axs[2].plot(df1.CZ[mskBASE_PH&(~mskDHI)],df1.DHI[mskBASE_PH&(~mskDHI)],'.b',label='pasan', markersize=2)
# axs[0].set(ylabel=r'$DHI(W/m^2)$',xlabel=r'$\cos(\theta_z)$',title=r'$DHI$',ylim=[-150,1500])
# axs[2].legend(markerscale=4, loc='best'); axs[2].grid()
# plt.tight_layout()


# fig, axs = plt.subplots(2, 3, figsize=(18,10))
# axs[0,0].plot(df1.CZ[mskBASE_PH&(mskGHI)],df1.GHI1[mskBASE_PH&(mskGHI)],'.',color='red', markersize=2)
# axs[0,0].set(ylabel='GHI(W/m^2)',xlabel='cos(theta_z)',title='GHI',ylim=[-100,1750])
# axs[0,0].grid()
# axs[0,1].plot(df1.CZ[mskBASE_PH&(mskDHI)],df1.DHI[mskBASE_PH&(mskDHI)],'.',color='red', markersize=2)
# axs[0,1].set(ylabel='DHI(W/m^2)',xlabel='cos(theta_z)',title='DHI',ylim=[-100,1750])
# axs[0,1].grid()
# axs[0,2].plot(df1.CZ[mskBASE_PH&(mskDNI)],df1.DNI1[mskBASE_PH&(mskDNI)],'.',color='red', markersize=2)
# axs[0,2].set(ylabel='DNI(W/m^2)',xlabel='cos(theta_z)',title='DNI',ylim=[-250,1750])
# axs[0,2].grid()
# axs[1,0].plot(df1.CZ[mskBASE_PH&(~mskGHI)],df1.GHI1[mskBASE_PH&(~mskGHI)],'.',color='blue', markersize=2)
# axs[1,0].set(ylabel='GHI(W/m^2)',xlabel='cos(theta_z)',title='GHI',ylim=[-100,1750])
# axs[1,0].grid()
# axs[1,1].plot(df1.CZ[mskBASE_PH&(~mskDHI)],df1.DHI[mskBASE_PH&(~mskDHI)],'.',color='blue', markersize=2)
# axs[1,1].set(ylabel='DHI(W/m^2)',xlabel='cos(theta_z)',title='DHI',ylim=[-100,1750])
# axs[1,1].grid()
# axs[1,2].plot(df1.CZ[mskBASE_PH&(~mskDNI)],df1.DNI1[mskBASE_PH&(~mskDNI)],'.',color='blue', markersize=2)
# axs[1,2].set(ylabel='DNI(W/m^2)',xlabel='cos(theta_z)',title='DNI',ylim=[-250,1750])
# axs[1,2].grid()
# #fig.suptitle("Datos que pasan el conjunto de filtros", fontsize=14)
# plt.tight_layout()

#%% REPORTE

flags = list(mskDict.keys())[1:]
reportePH = pd.DataFrame(index=flags,columns=['descarte','%desc'])

total = mskBASE_PH.sum() #total de datos válidos diurnos
totalPI = mskBASE_PI.sum() #total de datos válidos diurnos

for f in flags:
    reportePH.loc[f,'descarte'] = (mskDict[f]&mskBASE_PH).sum()
    reportePH.loc[f,'%desc'] = round(100 * ((mskDict[f]&mskBASE_PH).sum()) / total, 1)
    
print('Cantidad de datos válidos diurnos (PH): {}'.format(total)); print('')
print('Cantidad de datos válidos diurnos (PI): {}'.format(totalPI)); print('')
print(reportePH); print('')


subflags = [f for f in flags if f not in ['F13','F19','F20','F21','F22','F23']]
descartePH = ((mascarasResult[subflags].any(axis='columns'))&(mskBASE_PH)).sum()
print('Descarte final PH: {} datos; {:.1f}%'.format(descartePH,100*descartePH/total))
    
subflagsPI = [f for f in flags if f not in ['F13']]
descartePI = ((mascarasResult[subflagsPI].any(axis='columns'))&(mskBASE_PI)).sum()
print('Descarte final PI: {} datos; {:.1f}%'.format(descartePI,100*descartePI/total))
    


#%% Guardar salidas
try:
    os.mkdir(rutaOUT)
except: pass
try: 
    os.mkdir('{}{}/'.format(rutaOUT,est))
    
except: pass
try: 
    os.mkdir('{}{}/PX{:02d}/'.format(rutaOUT,est,PX))
except: pass


reportePH.to_csv('{}{}/PX{:02d}/PX{:02d}_{}_reportePH_v2.csv'.format(rutaOUT,est,PX,PX,est))

#%% Mostrar la proporción estacional de datos.

# Zthr = 80
# mskCZ = df1.CZ > np.cos(np.deg2rad(Zthr))
plt.rcParams.update({'font.size': 14})
plt.figure(); plt.title('fd-kt con mascaras de GHI y DHI ')
plt.plot(kt[mskBASE_PH],fd[mskBASE_PH],'.',color='grey', label='descarte', markersize=2)
plt.plot(kt[mskBASE_PH&(~mskGHI)&(~mskDHI)],fd[mskBASE_PH&(~mskGHI)&(~mskDHI)],'.k',label='pasan', markersize=2)
plt.xlabel(r'$k_t$'); plt.ylabel(r'$f_d$')
plt.xlim([-0.05,1.2]); plt.ylim([-0.05,1.2])
plt.legend(markerscale=4); plt.grid()
plt.tight_layout()
plt.savefig('{}/PX{:02d}_fd_kt_postfilt.png'.format(rutaFIG,PX))

# plt.figure(); plt.title('fd-kt / F16 mas z > 75° / TL = TLciclos - {}'.format(TL_01))
# plt.plot(kt[mskBASE_PH],fd[mskBASE_PH],'.',color='orange', label='Z > 75°', markersize=2)
# plt.plot(kt[mskBASE_PH&(mskF16)&mskCZ],fd[mskBASE_PH&(mskF16)&mskCZ],'.r',label='F16', markersize=2)
# plt.plot(kt[mskBASE_PH&(~mskF16)&mskCZ],fd[mskBASE_PH&(~mskF16)&mskCZ],'.b',label='pasan', markersize=2)
# plt.xlabel('kt'); plt.ylabel('fd')
# plt.xlim([-0.05,1.2]); plt.ylim([-0.05,1.2])
# plt.legend(markerscale=4,bbox_to_anchor=(1.04, 0.6)); plt.grid()
# plt.tight_layout()



# plt.figure(); plt.title('Descarte de filtro GTIcs')
# plt.plot(df1.index,GTIcs_0,'-r',label='GTIcs',linewidth=1)
# plt.plot(df1[mskBASE_PH].index,df1[mskBASE_PH].GTI,'.',color='grey',label='descarte')
# plt.plot(df1[mskBASE_PH&(~mskF19)&(~mskF20)].index,df1[mskBASE_PH&(~mskF19)&(~mskF20)].GTI,'.',color='black',label='pasan')
# plt.ylim([-100,2000])
# plt.legend(markerscale=4); plt.grid()

#%% EXTRAS
# mskGHI12 = abs((df1.GHI1-df1.GHI2)/df1.GHI1) > 0.1

# plt.figure(); plt.title(r'$f_d - k_t$ con mascaras de GHI y DHI ')
# plt.plot(kt[mskBASE_PH],fd[mskBASE_PH],'.',color='grey', label='descarte', markersize=2)
# plt.plot(kt[mskBASE_PH&(~mskGHI)&(~mskDHI)&(~mskGHI12)],fd[mskBASE_PH&(~mskGHI)&(~mskDHI)&(~mskGHI12)],'.k',label='pasan', markersize=2)
# plt.xlabel(r'$k_t$'); plt.ylabel(r'$f_d$')
# plt.xlim([-0.05,1.2]); plt.ylim([-0.05,1.2])
# plt.legend(markerscale=4); plt.grid()

# # #GHI1 vs GHI2
# plt.figure()
# plt.plot(df1.GHI1,df1.GHI2,'.',color='grey',markersize=2)
# plt.plot(df1[mskBASE_PH&(~mskGHI)&(~mskDHI)].GHI1,df1[mskBASE_PH&(~mskGHI)&(~mskDHI)].GHI2,'.k',markersize=2)
# plt.plot(df1[mskBASE_PH&(~mskGHI)&(~mskDHI)&(~mskGHI12)].GHI1,df1[mskBASE_PH&(~mskGHI)&(~mskDHI)&(~mskGHI12)].GHI2,'.g',markersize=2)
# plt.plot([0,1500],[0,1500],'-r')
# plt.xlim([-100,1500]); plt.ylim([-100, 1500])
# plt.xlabel(r'$GHI_1(W/m^2)$');plt.ylabel(r'$GHI_2(W/m^2)$')
# plt.grid()