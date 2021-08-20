#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 11:35:23 2021

@author: inti
"""
#Funciones para determinar el omega que maximiza la medida de GTI (ajustada)
import pandas as pd
import numpy as np
import pvlib as pv
from fun import indicadores
from fun import modelado as md

def getIntervalGTImax(GTI,deltaT=1): 
    #Devuelve la GTI alrededor del máximo (+- 1h por def)
    t_max = GTI.idxmax()
    t1 = t_max - pd.Timedelta(hours=deltaT)
    t2 = t_max + pd.Timedelta(hours=deltaT)
    return GTI.loc[t1:t2], t1, t2

def ajusteGTI(GTImax_intervalo,grado=4):
    #Devuelve la GTI ajustada por un pol. de gr 4 
    x_val = np.linspace(0,1,len(GTImax_intervalo))
    y_val = GTImax_intervalo
    
    coeffs = np.polyfit(x_val,y_val,4)
    poly_eqn = np.poly1d(coeffs)
    y_hat = pd.Series(poly_eqn(x_val),index=GTImax_intervalo.index)
    
    return y_hat

def getDeltaOmega(GTI_ajustada,omegaDeg):
    tmax = GTI_ajustada.idxmax()
    dw = omegaDeg[tmax]
    
    return dw

#Funciones para calcular la curva característica omega vs gamma.
def getTL(ciclosTL,est,doy):
    TL = ciclosTL[est]
    doy_medio = [15, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350] #doy del centro de cada mes
    TLinterp = np.interp(doy , doy_medio, TL, period = 366)
    
    return TLinterp
    

def getGTIcsk(BETAdeg,az,Zdeg,gammaDeg,Go,alturaSolar,Zsnm,TL=2.5):
    
    GHIesra,DHIesra,DNIesra = md.generaComponentesESRA(alturaSolar,Go,Zsnm,TL=TL)
    GTIesra = (pv.irradiance.get_total_irradiance(BETAdeg, az, Zdeg, gammaDeg, 
                                             DNIesra, GHIesra, DHIesra,Go, albedo=0.25))
    GTIesra = GTIesra['poa_global']
    
    return GTIesra


def getSerieTempAux(datos):
    #Creo una serie temporal con resulución de segundos, para el día de esos datos
    serieTemp_1s = datos.asfreq('1s').index
    return serieTemp_1s

def calculosSolares(loc,serieTemp_1s):
    #Calculo el df solar_pos que guarda información de los ángulos solares
    solar_pos_1s = loc.get_solarposition(serieTemp_1s)
    
    return solar_pos_1s

def curvaCaracteristica(loc,solar_pos_1s,BETAdeg,Zsnm,TL,azimuts):
    ''' Calcula las delta omega para el rango dado por (-az_lim,az_lim). Devuelve
    el rango de azimuts y las delta omega correspondientes
    '''
    #Extraigo/calculo las variables solares relevantes:
    Zdeg = solar_pos_1s['zenith']
    gammaDeg = solar_pos_1s['azimuth']
    alturaSolar = solar_pos_1s['elevation']
    omegaDeg = pv.solarposition.hour_angle(solar_pos_1s.index, loc.longitude, solar_pos_1s.equation_of_time)
    omegaDeg = pd.Series(omegaDeg,index = solar_pos_1s.index)
    Go = pv.irradiance.get_extra_radiation(solar_pos_1s.index)
  
    #Recorre el rango de azimuts, calcula GTIesra y determina el omega que maximiza GTIesra
    dw_azimutal = np.zeros(len(azimuts))
    j = 0
    for az_deg in azimuts:
        GTIaux = getGTIcsk(BETAdeg,az_deg,Zdeg,gammaDeg,Go,alturaSolar,Zsnm,TL)
        dw_azimutal[j] = omegaDeg.loc[GTIaux.idxmax()]
        j += 1
        
    
    return dw_azimutal

def ajustaCurva(azimuts,dw_azimutal):
    '''Ajusta la curva delta omega en función del azimut por un polinomio de 3º grado, devolviendo los valores
    funcionales del ajuste así como la ecuación del mismo.
    '''
    x_val = azimuts; y_val = dw_azimutal
    coeffs = np.polyfit(x_val,y_val,3)
    poly_eqn = np.poly1d(coeffs)
    y_hat = poly_eqn(x_val)
    
    return y_hat, poly_eqn

def interpola(dwExperimental,poly_eqn,azimuts):
    '''Recibe un delta omega obtenido experimentalmente y el resultado del ajuste de la curva, interpola
    y devuelve el resultado, que corresponde a la estimación del azimut (dado que es el corte de una recta
    con un polinomio de 3º grado, son 3 soluciones y hay que elegir la que está en el intervalo (-az,az))
    '''
    raices = np.roots(poly_eqn - dwExperimental)
    raiz = raices[(raices<azimuts[-1])&(raices>(azimuts[0]))]
    
    return raiz

#Funciones misceláneas:
def preparaDatos(rutaPX,year,PX,est,doy,tz):
    if (est=='TT') | (est=='TA'):
        fileName = rutaPX + '{0}/PX{1:0=2d}_{2}_RAD_{3}{4:0=3d}.csv'.format(year,PX,est,year,doy)
    elif est=='AR':
        fileName = rutaPX + '{0}/PX{1:0=2d}_{2}c_{3}{4:0=3d}.csv'.format(year,PX,est,year,doy)
    datos = pd.read_csv(fileName, index_col = 'Fecha')
    datos.index = pd.to_datetime(datos.index)
    datos.index = datos.index.tz_localize(tz)
    
    datos = datos[['GHI1','GTI']] #me quedo  con las columnas relevantes
    
    return datos

def getVariablesSolares(solar_pos,loc):
    Zdeg = solar_pos['zenith']
    gammaDeg = solar_pos['azimuth']
    alturaSolar = solar_pos['elevation']
    omegaDeg = pv.solarposition.hour_angle(solar_pos.index, loc.longitude, solar_pos.equation_of_time)
    omegaDeg = pd.Series(omegaDeg,index = solar_pos.index)
    Go = pv.irradiance.get_extra_radiation(solar_pos.index)
    
    return Zdeg, gammaDeg, alturaSolar, omegaDeg, Go

#%% Funciones para el segundo método de detección de azimut (mínimo de indicadores del pasaje a PI con azimut variable)

def getAZmetodo2(BETAdeg,zenith,azimutSolar,DNI,GHI,DHI,GTI,AM,Go,azimuts,mskFINAL,albedo=0.2, model='isotropic'):
    rms = np.zeros(len(azimuts))
    i = 0
    for gamma in azimuts:
        GTIm = (pv.irradiance.get_total_irradiance(BETAdeg,gamma,zenith[mskFINAL],azimutSolar[mskFINAL],
                                        DNI[mskFINAL],GHI[mskFINAL],DHI[mskFINAL],Go[mskFINAL],AM[mskFINAL],albedo=albedo,model=model))
        GTIm = GTIm['poa_global']
        GTImedia,_,_,rMBE,rRMSD = indicadores.desvios(GTIm[mskFINAL], GTI[mskFINAL])
            
        
        rms[i] = rRMSD 
        
        i = i + 1
    
    AZ = azimuts[rms==rms.min()][0]
    
    return AZ, rms