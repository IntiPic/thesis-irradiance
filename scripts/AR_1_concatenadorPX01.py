#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 09:48:40 2021

CONCATENADOR PX01 

@author: inti
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pvlib as pv
import datetime as dt
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
#%% Creación de un df vacío con el total de los datos
years = range(2017,2019)
t_ini = pd.Timestamp('1/1/{}'.format(years[0])); t_fin = pd.Timestamp('12/31/{} 23:59:00'.format(years[-1]))
serieTemp = pd.date_range(start=t_ini, end=t_fin, freq="{:02d}".format(PX) + 'min', tz=tz)

cols =  ['CZ','GHI1', 'GHI2', 'DHI', 'GTI']

dfConcat = pd.DataFrame(np.nan, index = serieTemp, columns = cols, dtype='float64')
dfConcat.index.name = 'Fecha'



for year in years:
    folderYear = rutaPX + '/{}/CSV/'.format(year)
    filesYear = sorted(os.listdir(folderYear))
    for file in filesYear:
        df01Year = pd.read_csv(folderYear + file, index_col='Fecha')
        df01Year.index = pd.to_datetime(df01Year.index)
        df01Year.index = df01Year.index.tz_localize(tz)
        

        dfConcat.loc[df01Year.index,cols[1:]] = df01Year
        print(file)

# try:
#     os.mkdir(rutaPX + '/concat')      #NO QUIERO GUARDAR EN UNA CARPETA "CONCAT" - REVISAR!
# except: pass

solar_pos = loc.get_solarposition(dfConcat.index)
Zdeg = solar_pos['zenith']

dfConcat['CZ'] = np.cos(np.deg2rad(Zdeg))

#Guardo la salida
dfConcat.to_csv('{}/PX{:02d}_{}_{}_{}.csv'.format(rutaPX,PX,est,years[0],years[-1]))

#%% PRUEBA GRÁFICA: serie temporal de cada columna:
    


plt.rcParams.update({'font.size': 16})
fig, axs = plt.subplots(4, 1, figsize=(18,18))
i=0

for col in cols[1:]:
    axs[i].plot(dfConcat[col])
    axs[i].set(ylabel='{}(W/m^2)'.format(col),xlabel='Fecha',ylim=[-100,1750])
    axs[i].grid()    
    i = i+1
    
#Guardo salida
#plt.savefig(rutaFIG+'/PX{:02d}_serieMinutal.png'.format(PX))
