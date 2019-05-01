import math
from .tempProfile import tempProfile
from .profile import *
class salinityProfile(Profile):
    def __init__(self,pressures,temperatures,salinities,densities,tp=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        if tp == None:
            tp = tempProfile(pressures,temperatures)
        tp.findMLD()
        self.MLDT = int(tp.foundMLD)
        self.dT = tp.dT
        ##fnd reference pressure, this is done in holte and talley 
        ##supplementary matlab file
        startindex = np.argmin((np.asarray(pressures)-10)**2)
        self.salinities = np.round_(salinities[startindex:],4)
        self.pressures = pressures[startindex:]
        self.densities = densities[startindex:]
        self.densityGradients = self.generateGradientList(self.densities)
        self.salinityGradients = np.round_(self.generateGradientList(self.salinities),5)
        self.SGradientMax = int(self.calculateSGradientMax())
        self.SGradientMaxPressure = self.pressures[self.SGradientMax]
        self.MLTFITSalinity, self.MLTFITSalinityPressure = self.calculateMLTFIT(self.salinities,self.salinityGradients)
        self.MLTFITDensity, self.MLTFITDensityPressure = self.calculateMLTFIT(self.densities,self.densityGradients)
        self.DThreshold = int(self.calculateDThreshold())
        self.DThresholdPressure = math.floor(self.interpolateDThreshold()*2)/2.0
        self.densityTest = self.calculateDensityTest()
        self.intrusionDepth = int(self.calculateIntrusionDepth())
        self.intrusionDepthPressure = self.pressures[self.intrusionDepth]
        self.foundMLD=0
        self.path=""
        self.debug=""
        #range from paper
        self.range = 25
        return

       #The minimum of the depth of the salinity gradiet maximum and the salinity minimum
    #In matlab file: dsandsmin
    def calculateIntrusionDepth(self):
        x = np.argwhere(self.salinityGradients == np.min(self.salinityGradients))[-1][0]+1
        y =  np.argwhere(self.salinities == np.min(self.salinities))[-1][0]
        if abs(self.pressures[y] - self.pressures[x] ) < 100:
            return min(y,x)
        else:
            return 0

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldWinterProfile(self):
        if self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
            self.debug="intrusionDepth zip"
            if MLD > self.DThresholdPressure:
                MLD = self.DThresholdPressure
                self.debug="DThreshold zip"

        else:
            if self.SGradientMaxPressure < self.DThresholdPressure:
                MLD = self.SGradientMaxPressure
                self.debug="SGradientMAx zip"
                if self.MLTFITSalinityPressure < MLD:
                    MLD = self.MLTFITSalinityPressure
                    self.debug="MLTFIT zip"
            else:
                MLD = self.DThresholdPressure
                self.debug = "DThreshold zap"
                if self.MLTFITSalinityPressure < MLD:
                    MLD = self.MLTFITSalinityPressure
                    self.debug = "MLTFIT zap"
                if MLD ==0:
                    MLD = self.SGradientMaxPressure
                    self.debug = "SGRadientMax zap"
                if self.SGradientMaxPressure > self.DThresholdPressure:
                    MLD = self.DThresholdPressure
                    self.debug = "DThreshold zop"
        return MLD

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldSummerProfile(self):
        MLD = self.MLTFITSalinityPressure
        self.debug = "MLTFIT zoop"
        if MLD - self.DThresholdPressure > self.range:
            MLD = self.DThresholdPressure
            self.debug = "DThreshold zoop"
        if self.MLTFITSalinityPressure - self.SGradientMaxPressure < 0 and self.DThresholdPressure - self.SGradientMaxPressure > 0:
            MLD = self.SGradientMaxPressure
            self.debug = "SGRadientMax zop"
        if self.MLTFITSalinityPressure - self.intrusionDepthPressure < self.range and self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
            self.debug = "SGRadientMax zeepp"
        if abs(self.DThresholdPressure - self.intrusionDepthPressure) < self.range and self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
            self.debug = "IntrusionDeth zinc"
        if self.MLDT - self.DThresholdPressure < 0 and abs(self.MLDT - self.DThresholdPressure) < self.range:
            MLD = self.MLDT
            self.debug = "MLDT"
            if abs(self.MLDT-self.MLTFITSalinityPressure)<self.range and self.MLTFITSalinityPressure - self.DThresholdPressure < 0:
                MLD = self.MLTFITSalinityPressure
                self.debug = "MLTFIT zop"
        if abs(self.MLDT - self.DThresholdPressure) < abs(MLD - self.DThresholdPressure):
            if self.MLDT > self.DThresholdPressure:
                MLD = self.DThresholdPressure
                self.debug = "DThreshold zoot"
        return MLD

    def findMLD(self):
        if self.densityTest == 1:
            self.foundMLD =  self.mldWinterProfile()
        else:
            self.foundMLD = self.mldSummerProfile()
        return self.foundMLD

    def importantDepths(self):
        return [
            [self.MLTFITSalinity,"MLTFIT"],
            [self.DThreshold,"Density Threshold"],
            [self.SGradientMax,"Salinity Gradient Maximum"],
            [self.intrusionDepth,"Intrusion Depth"]
        ]

    def __str__(self):
        out = ""
        out += "MLTFITSalinityPressure: " + str(self.MLTFITPressure) + "\n"
        out += "DThresholdPressure: " + str(self.DThresholdPressure) + "\n"
        out += "SGradientMaxPressure: " + str(self.SGradientMaxPressure) + "\n"
        out += "intrusionDepthPressure: " + str(self.intrusionDepthPressure) + "\n"
        out += "MLD Depth: " + str(self.foundMLD) + "\n"
        return out
