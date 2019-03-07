import math
import numpy as np
from .profile import *
from .tempProfile import tempProfile
from .salinityProfile import salinityProfile
class densityProfile(Profile):
    def __init__(self,pressures,temperatures,salinities,densities,tp=None,sp=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        if tp == None:
            self.tp = tempProfile(pressures,temperatures)
        else:
            self.tp = tp
        if sp == None:
            self.sp = salinityProfile(pressures,temperatures,salinities,densities)
        else:
            self.sp = sp
        self.tp.findMLD()
        self.MLDT = self.tp.foundMLD
        self.sp.findMLD()
        self.MLDS = self.sp.foundMLD
        ##fnd reference pressure, this is done in holte and talley 
        ##supplementary matlab file
        startindex = np.argmin((np.asarray(pressures)-10.0)**2)
        self.salinities = salinities[startindex:]
        self.pressures = pressures[startindex:]
        self.densities = densities[startindex:]
        self.densityGradients = self.generateGradientList(self.densities)
        self.DMin = int(self.calculateDMin())
        self.DMinPressure = pressures[self.DMin]
        self.MLTFITDensity, self.MLTFITDensityPressure  = self.calculateMLTFIT(self.densities,self.densityGradients)
        self.DThreshold = int(self.calculateDThreshold())
        self.DThresholdPressure = pressures[self.DThreshold]
        self.DThresholdPressure = math.floor(self.interpolateDThreshold()*2)/2.0
        self.DGradientThreshold = int(self.calculateDGradientThreshold())
        self.DGradientThresholdPressure = pressures[self.DGradientThreshold]
        self.densityTest = sp.calculateDensityTest()
        self.foundMLD=0
        self.path=""
        self.debug=""
        #range from paper
        self.range = 25
        return

    #The Density Minimum
    #Closest paper equivalent TMAX
    def calculateDMin(self):
        return np.argmin(self.densities)

    #Calculates the densityGradient threshold or max if criteria not met
    #closest to DTM
    def calculateDGradientThreshold(self):
        for index in range(0,len(self.densityGradients)):
            if abs(self.densityGradients[index]) > 0.0005:
                return index+2
        return np.argmax(self.densityGradients)+1

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldWinterProfile(self):
        MLD = self.DThresholdPressure
        self.debug = "D threshold zip"
        if self.tp.TTMLDPressure < MLD:
            MLD = self.tp.TTMLDPressure
            self.debug = "TTMLD zap"
        if self.MLTFITDensityPressure < self.DThresholdPressure and self.MLTFITDensityPressure > self.range:
            MLD = self.MLTFITDensityPressure
            self.debug = "MLTFIT zip"
        if self.tp.TDTMPressure > self.range and self.tp.TDTMPressure < self.DThresholdPressure:
            MLD= self.tp.TDTMPressure
            self.debug = "TDTM zip"
            if abs(self.tp.TMaxPressure - self.DThresholdPressure) < abs(self.tp.TDTMPressure - self.MLTFITDensityPressure):
                MLD = self.tp.TMaxPressure
                self.debug = "TMAX zip"
            if abs(self.MLDS - self.DThresholdPressure) < self.range and self.MLDS < self.DThresholdPressure:
                MLD = min(self.DThresholdPressure,self.MLDS)
                self.debug = "DThreshold zap"
        if abs(self.MLDT - self.MLDS) < self.range and abs(min(self.MLDT,self.MLDS))-MLD > self.range:
            MLD = min(self.MLDT,self.MLDS)
            self.debug = "Min of MLDT and MLDS zip"
        if MLD > self.DGradientThresholdPressure and abs(self.DGradientThresholdPressure - self.MLDT) < abs(MLD-self.MLDT):
            MLD = self.DGradientThresholdPressure
            self.debug = "D Gradient Threshold zip"
        if self.MLTFITDensityPressure == self.sp.MLTFITSalinityPressure and abs(self.sp.MLTFITSalinityPressure - self.DThresholdPressure) < self.range:
            MLD= self.MLTFITDensityPressure
            self.debug = "MLTIT zap"
        if self.MLDT == self.DMinPressure:
            MLD = self.DMinPressure
            self.debug = "Dmin zip"
        return MLD

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldSummerProfile(self):
        MLD = self.MLTFITDensityPressure
        self.debug = "MLTFIT zip"
        if MLD > self.DThresholdPressure:
            MLD = self.DThresholdPressure
            self.debug = "D treshold zap"
        ## If any two of these three conditions are true
        conditions = [abs(self.MLDS-self.MLDT),abs(self.MLTFITDensityPressure-self.MLDT),abs(self.MLDS - self.MLTFITDensityPressure)]
        s =0
        for a in conditions:
            if a <self.range:
                s+=1
        if s >1:
            MLD = self.MLTFITDensityPressure
            self.debug = "MLTFIT zap"
        if abs(self.MLDS - self.DThresholdPressure) < self.range and self.MLDS != self.DThresholdPressure:
            if self.DThresholdPressure < self.MLDS:
                MLD = self.DThresholdPressure
                self.debug = "D treshold zop"
            else:
                MLD = self.MLDS
                self.debug = "MLDS zip"
            if self.MLTFITDensityPressure == self.DThresholdPressure:
                MLD = self.MLTFITDensityPressure
                self.debug = "MLTFIT zap"
        if MLD > self.DGradientThresholdPressure and abs(self.DGradientThresholdPressure - self.MLDT) < abs(MLD - self.MLDT):
            MLD = self.DGradientThresholdPressure 
            self.debug = "D Gradient Threshold zap"
        return MLD

    def findMLD(self):
        if self.sp.densityTest:
            self.foundMLD =  self.mldWinterProfile()
        else:
            self.foundMLD = self.mldSummerProfile()
        return self.foundMLD

    def importantDepths(self):
        return [
            [int(self.MLTFIT),"mltfit pressure"],
            #[],
            #[],
            #[],
            #[]
        ]

    def __str__(self):
        out = ""
        out += "MLDD: " + str(self.foundMLD) + "\n"
        out += "MLDS: " + str(self.MLDS) + "\n"
        out += "MLDT: " + str(self.MLDT) + "\n"
        out += "" + str() + "\n"
        out += "" + str() + "\n"
        out += "" + str() + "\n"
        return out
