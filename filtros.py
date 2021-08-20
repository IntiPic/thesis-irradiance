#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 11:50:20 2021

@author: inti
"""

import numpy as np
import pandas as pd


def F0(GHI,GHImin=-4):
    return GHI < GHImin

def F1(DHI,DHImin=-4):
    return DHI < DHImin

def F2(DNI,DNImin=-4):
    return DNI < DNImin

def F3(GHI,CZ,Iob,a=1.5,b=1.2,c=100):
    return GHI > Iob*a*(CZ**b) + c
    
def F4(DHI,CZ,Iob,a=0.95,b=1.2,c=50):
    return DHI > Iob*a*(CZ**b) + c
    
def F5(DNI,Iob):
    return DNI > Iob

def F6(DNI,DNIcs): #DNIcs: DNI de cielo claro, según modelo ESRA, TL=1
    return DNI > DNIcs

def F7(GHI,GHIsum,CZ,delta=0.08,GHIsum_tol = 50):
    theta_z = np.rad2deg(np.arccos(CZ))
    return (abs(GHI/GHIsum - 1) > delta) & (theta_z <= 75) & (GHIsum > GHIsum_tol)
    
def F8(GHI,GHIsum,CZ,delta=0.15,GHIsum_tol = 50):
    theta_z = np.rad2deg(np.arccos(CZ))
    return (abs(GHI/GHIsum - 1) > delta) & ((theta_z > 75)&(theta_z < 93)) & (GHIsum > GHIsum_tol)

def F9(GHI,DHI,CZ,fdmax = 1.05,GHImin=50):
    fd = DHI/GHI
    theta_z = np.rad2deg(np.arccos(CZ))
    return (fd > fdmax) & (theta_z <= 75) & (GHI > GHImin)

def F10(GHI,DHI,CZ,fdmax = 1.10,GHImin=50):
    fd = DHI/GHI
    theta_z = np.rad2deg(np.arccos(CZ))
    return (fd > fdmax) & ((theta_z > 75)&(theta_z < 93)) & (GHI > GHImin)

def F11(GHI,DHI,DNI,CZ,tol=50):
    return abs(DNI*CZ - (GHI - DHI)) > tol

def F12(DHI,DHImax=700):
    return DHI > DHImax

def F13(DHI,Iob,CZ,tol=0.6):
    return DHI/(Iob*CZ) > tol

def F14(GHI,DHI,CZ,Iob,ktmax=0.2,fdmax=0.9):
    fd = DHI/GHI
    kt = GHI/(Iob*CZ)
    return (kt < ktmax) & (fd<fdmax)

def F15(GHI,DHI,CZ,Iob,ktmin=0.5,fdmin=0.8):
    fd = DHI/GHI
    kt = GHI/(Iob*CZ)
    return (kt > ktmin) & (fd > fdmin)

def F16(DHI,GHIsum,GHIcs,DHImin=50,tol=0.85):
    return (GHIsum/GHIcs > tol) & (DHI/GHIsum > tol) & (DHI > DHImin)
    
def F17(GHI,Iob,CZ,AM,tol=1):
    kt = GHI/(Iob*CZ)
    kt_per = kt/(1.031*np.exp(-1.4/(0.9 + (9.4/AM))) + 0.1)
    return kt_per > tol

def F18(GHI,DHI,DNIcs,CZ):
    return ((GHI - DHI)/CZ) > DNIcs

    

#%% EXTRA: balance de datos estacional
def proporcionEstacional(mskTOTAL):
    serieTemp = mskTOTAL.index
    mskVER = (serieTemp.dayofyear > 355)|(serieTemp.dayofyear <= 80)
    mskINV = (serieTemp.dayofyear > 172)&(serieTemp.dayofyear <= 264)
    mskME = (((serieTemp.dayofyear > 264)&(serieTemp.dayofyear <= 355))|
            (serieTemp.dayofyear > 80)&(serieTemp.dayofyear <= 172))
    
    df = pd.Series(index = ['%VER','%INV','%ME'])
    df['%VER'] = 100*(mskVER&mskTOTAL).sum()/mskTOTAL.sum()
    df['%INV'] = 100*(mskINV&mskTOTAL).sum()/mskTOTAL.sum()
    df['%ME'] = 100*(mskME&mskTOTAL).sum()/mskTOTAL.sum()
     
    return df.astype(float).round(1)


