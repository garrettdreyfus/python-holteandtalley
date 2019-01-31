import numpy as np
import math
from tempProfile import tempProfile
class salinityProfile:
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
        self.steepest = np.argwhere(np.abs(self.salinityGradients) == np.max(np.abs(self.salinityGradients)))[-1][0]+1
        self.MLTFIT, self.MLTFITPressure = self.calculateMLTFIT(self.salinities,self.salinityGradients)
        self.MLTFITDensity, self.MLTFITDensityPressure = self.calculateMLTFIT(self.densities,self.densityGradients)
        self.DThreshold = int(self.calculateDThreshold())
        self.DThresholdPressure = math.floor(self.interpolateDThreshold()*2)/2.0
        self.densityTest = self.calculateDensityTest()
        self.SGradientMax = int(self.calculateSGradientMax())
        self.SGradientMaxPressure = self.pressures[self.SGradientMax]
        self.intrusionDepth = int(self.calculateIntrusionDepth())
        self.intrusionDepthPressure = self.pressures[self.intrusionDepth]
        self.foundMLD=0
        self.path=""
        self.debug=""
        #range from paper
        self.range = 25
        return

    ##find nearest pressure
    def findNearestPressureIndex(self,value):
        return (np.abs(np.asarray(self.pressures)- value)).argmin()

    def generateGradientList(self,values):
        tGS=[]
        for index in range(len(values)-1):
            dt = float(values[index] - values[index+1])
            dp = float(self.pressures[index] - self.pressures[index+1])
            tGS.append(dt/dp)
        smoothed=[0]*(len(tGS)-2)
        for i in range(1,len(tGS)-1):
            smoothed[i-1] = (tGS[i-1]+tGS[i]+tGS[i+1])/3.0
        return smoothed

    #calculates the TTMLD or the salinity threshold mixed layer estimate
    # based on a salinity threshold of 0.03 from Boyer Montegut
    # Matlab equivalent: mldepthdens
    #closest paper equivalent: TTMLD
    def calculateDThreshold(self):
        for index in range(0,len(self.pressures)):
            if abs(self.densities[index] - self.densities[0]) > 0.03:
                return index
        return 0
    def interpolateDThreshold(self):
        if self.DThreshold:
            thresholdIndex = self.DThreshold
        else:
            thresholdIndex = int(self.calculateDThreshold())
        deltaP = (self.pressures[thresholdIndex] - self.pressures[thresholdIndex-1])
        deltaD = (self.densities[thresholdIndex] - self.densities[thresholdIndex-1])
        if self.densities[thresholdIndex] > self.densities[thresholdIndex-1]:
            return self.pressures[thresholdIndex-1] + (deltaP/deltaD)*(self.densities[0]+0.03-abs(self.densities[thresholdIndex-1]))
        elif self.densities[thresholdIndex] < self.densities[thresholdIndex-1]:
            return self.pressures[thresholdIndex-1] + (deltaP/deltaD)*(self.densities[0]-0.03-abs(self.densities[thresholdIndex-1]))
        #return self.pressures[thresholdIndex]
    # calculates the MLTFIT or the intersection of the best fit mixed layer line
    # and the best fit of the thermocline
    #Matlab equivalent: upperdsmax
    #closest paper equivalent: MLTFIT
    def calculateMLTFIT(self,values,gradients):
        #Calculate the best fit line of the mixed layer
        errors = []
        fits = []
        #iterate through and polyfit over progressively increasing points
        for num in range(2,len(self.pressures)):
            out = np.polyfit(self.pressures[0:num],values[0:num],1,full=True)
            fits.append(out[0])
            if out[1]:
                errors.append(out[1][0])
            else:
                errors.append(0)
        errorsum = np.sum(errors)
        #find line with normalized error less than 10^-10
        self.errors= errors
        mltBestFit=None
        mltBestFit = np.max(np.argwhere((errors/errorsum)> 10**-10))
        for index in range(len(errors)):
            if errors[index]/errorsum >(10**-10):
                mltBestFit = fits[index-1]
                break
        self.mltfitline=mltBestFit
        self.mltfitindex = index-1
        #find thermoclineFit
        steepest = np.argwhere(np.abs(gradients) == np.max(np.abs(gradients)))[-1][0] +1
        #thermoclineFit = [gradients[steepest],
            #values[steepest]-gradients[steepest]*self.pressures[steepest]
        #]
        thermoclineFit = np.polyfit(self.pressures[steepest-1:steepest+2],values[steepest-1:steepest+2],1,full=True)[0]
        depth = abs(float(thermoclineFit[1] - mltBestFit[1])/float(thermoclineFit[0] - mltBestFit[0]))
        if False:
            return self.findNearestPressureIndex(depth), depth
        else:
            return self.findNearestPressureIndex(depth),self.pressures[self.findNearestPressureIndex(depth)]# depth

    # TESTD from matlab file
    def calculateDensityTest(self):
        if self.MLTFITDensity > 0 and self.MLTFITDensity < len(self.pressures)-2:
            ddiff = self.densities[self.MLTFITDensity]-self.densities[self.MLTFITDensity+2]
        else:
            densityGradientMax = np.argmax(np.abs(self.densityGradients))+1
            ddiff = self.densities[densityGradientMax-1] - self.densities[densityGradientMax+1]
        #various constants from paper
        if ddiff > -0.06 and self.dT > 0.5:
            return 1
        if ddiff > -0.06 and self.dT < -0.25:
            return 0
        if self.dT > 0.5 or self.dT < -0.25:
            self.ddiff = [ddiff, self.dT]
            return 0 
        else:
            return 1

    #The salinity gradient threshold or max if threshold not met
    # Matlab: dsmin (confusingly named I know)
    def calculateSGradientMax(self):
        return self.steepest

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
                if self.MLTFITPressure < MLD:
                    MLD = self.MLTFITPressure
                    self.debug="MLTFIT zip"
            else:
                MLD = self.DThresholdPressure
                self.debug = "DThreshold zap"
                if self.MLTFITPressure < MLD:
                    MLD = self.MLTFITPressure
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
        MLD = self.MLTFITPressure
        self.debug = "MLTFIT zoop"
        if MLD - self.DThresholdPressure > self.range:
            MLD = self.DThresholdPressure
            self.debug = "DThreshold zoop"
        if self.MLTFITPressure - self.SGradientMaxPressure < 0 and self.DThresholdPressure - self.SGradientMaxPressure > 0:
            MLD = self.SGradientMaxPressure
            self.debug = "SGRadientMax zop"
        if self.MLTFITPressure - self.intrusionDepthPressure < self.range and self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
            self.debug = "SGRadientMax zeepp"
        if abs(self.DThresholdPressure - self.intrusionDepthPressure) < self.range and self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
            self.debug = "IntrusionDeth zinc"
        if self.MLDT - self.DThresholdPressure < 0 and abs(self.MLDT - self.DThresholdPressure) < self.range:
            MLD = self.MLDT
            self.debug = "MLDT"
            if abs(self.MLDT-self.MLTFITPressure)<self.range and self.MLTFITPressure - self.DThresholdPressure < 0:
                MLD = self.MLTFITPressure
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
            [self.MLTFIT,"MLTFIT"],
            [self.DThreshold,"Density Threshold"],
            [self.SGradientMax,"Salinity Gradient Maximum"],
            [self.intrusionDepth,"Intrusion Depth"]
        ]

    def __str__(self):
        out = ""
        out += "MLTFITPressure: " + str(self.MLTFITPressure) + "\n"
        out += "DThresholdPressure: " + str(self.DThresholdPressure) + "\n"
        out += "SGradientMaxPressure: " + str(self.SGradientMaxPressure) + "\n"
        out += "intrusionDepthPressure: " + str(self.intrusionDepthPressure) + "\n"
        out += "MLD Depth: " + str(self.foundMLD) + "\n"
        return out
