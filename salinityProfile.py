import numpy as np
from tempProfile import tempProfile
class salinityProfile:
    def __init__(self,pressures,temperatures,salinities,densities,tp=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        if tp == None:
            tp = tempProfile(pressures,temperatures)
        tp.findMLD()
        self.MLDT = tp.foundMLD
        self.dT = tp.dT
        ##fnd reference pressure, this is done in holte and talley 
        ##supplementary matlab file
        startindex = np.argmin((np.asarray(pressures)-10)**2)
        self.salinities = salinities[startindex:]
        self.pressures = pressures[startindex:]
        self.densities = densities[startindex:]
        self.salinityGradients = self.generateGradientList(self.salinities)
        self.densityGradients = self.generateGradientList(self.densities)
        self.SMin = int(self.calculateSMin())
        self.SMinPressure = pressures[self.SMin]
        self.MLTFIT = int(self.calculateMLTFIT(self.salinities,self.salinityGradients))
        self.MLTFITDensity = int(self.calculateMLTFIT(self.densities,self.densityGradients))
        self.MLTFITPressure = pressures[self.MLTFIT]
        self.DThreshold = int(self.calculateDThreshold())
        self.DThresholdPressure = pressures[self.DThreshold]
        self.densityTest = self.calculateDensityTest()
        self.SGradientMax = int(self.calculateSGradientMax())
        self.SGradientMaxPressure = self.pressures[self.SGradientMax]
        self.intrusionDepth = int(self.calculateIntrusionDepth())
        self.intrusionDepthPressure = self.pressures[self.intrusionDepth]
        self.foundMLD=0
        self.path=""
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

    #The salinity minimum
    #Closest paper equivalent TMAX
    def calculateSMin(self):
        maxIndex = 0
        for i in range(len(self.salinities)):
            if self.salinities[i] >= self.salinities[maxIndex]:
                maxIndex = i
        return maxIndex

    #calculates the TTMLD or the salinity threshold mixed layer estimate
    # based on a salinity threshold of 0.03 from Boyer Montegut
    # Matlab equivalent: mldepthdens
    #closest paper equivalent: TTMLD
    def calculateDThreshold(self):
        for index in range(0,len(self.pressures)):
            if abs(self.densities[index] - self.densities[0]) > 0.03:
                return index
        return 0

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
        mltBestFit=None
        for index in range(len(errors)):
            if errors[index]/errorsum >(10**-10):
               mltBestFit = fits[index-1]
               break
        self.mltfitline=mltBestFit
        #find thermoclineFit
        steepest = np.argmax(np.abs(gradients))
        thermoclineFit = [gradients[steepest],
            values[steepest]-gradients[steepest]*self.pressures[steepest]
        ]
        self.thermoclinefitline = thermoclineFit
        depth = abs(float(thermoclineFit[1] - mltBestFit[1])/float(thermoclineFit[0] - mltBestFit[0]))
        return self.findNearestPressureIndex(depth)

    # TESTD from matlab file
    def calculateDensityTest(self):
        if self.MLTFITDensity > 0 and self.MLTFITDensity < len(self.pressures)-2:
            ddiff = self.densities[self.MLTFIT]-self.densities[self.MLTFIT+2]
        else:
            densityGradientMax = np.argmax(self.densityGradients)
            ddiff = self.densities[densityGradientMax-1] - self.densities[densityGradientMax+1]
        #various constants from paper
        if ddiff > -0.06 and self.dT > 0.5:
            return 1
        if ddiff > -0.06 and self.dT < -0.25:
            return 0
        if self.dT > -0.25 and self.dT < 0.5:
            return 1 
        else:
            return 0

    #The salinity gradient threshold or max if threshold not met
    # Matlab: dsmin (confusingly named I know)
    def calculateSGradientMax(self):
        return np.argmax(self.salinityGradients)+1

    #The minimum of the depth of the salinity gradiet maximum and the salinity minimum
    #In matlab file: dsandsmin
    def calculateIntrusionDepth(self):
        steepest=0
        for i in range(len(self.salinityGradients)):
            if  (self.salinityGradients[i]) > self.salinityGradients[steepest]:
                steepest=i
        if self.pressures[steepest] < self.pressures[self.SMin]:
            return steepest
        else:
            return self.SMin

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldWinterProfile(self):
        if self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
            if MLD > self.DThresholdPressure:
                MLD = self.DThresholdPressure
        else:
            if self.SGradientMaxPressure < self.DThresholdPressure:
                MLD = self.SGradientMaxPressure
                if self.MLTFITPressure < MLD:
                    MLD = self.MLTFITPressure
            else:
                MLD = self.DThresholdPressure
                if self.MLTFITPressure < MLD:
                    MLD = self.MLTFITPressure
                if MLD ==0:
                    MLD = self.SGradientMaxPressure
                if self.SGradientMaxPressure > self.DThresholdPressure:
                    MLD = self.DThresholdPressure
        return MLD

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldSummerProfile(self):
        MLD = self.MLTFIT
        if MLD - self.DThresholdPressure > self.range:
            MLD = self.DThresholdPressure
        if self.MLTFITPressure - self.SGradientMaxPressure < 0 and self.DThresholdPressure - self.SGradientMaxPressure > 0:
            MLD = self.SGradientMaxPressure
        if abs(self.MLTFITPressure - self.intrusionDepthPressure) < self.range and self.intrusionDepthPressure > self.range:
            MLD = self.intrusionDepthPressure
        if self.MLDT - self.DThresholdPressure < 0 and abs(self.MLDT - self.DThresholdPressure) < self.range:
            MLD = self.MLDT
            if abs(self.MLDT-self.MLTFITPressure)<self.range and self.MLTFITPressure - self.DThresholdPressure < 0:
                MLD = self.MLTFITPressure
        if abs(self.MLDT - self.DThresholdPressure) < abs(MLD - self.DThresholdPressure):
            if self.MLDT > self.DThresholdPressure:
                MLD = self.DThresholdPressure
        return MLD

    def findMLD(self):
        if self.densityTest:
            self.foundMLD =  self.mldWinterProfile()
        else:
            self.foundMLD = self.mldSummerProfile()
        return self.foundMLD

    def importantDepths(self):
        return [
            [self.SMin,"Salinity Minimum"],
            [self.MLTFIT,"MLTFIT"],
            [self.DThreshold,"Density Threshold"],
            [self.SGradientMax,"Salinity Gradient Maximum"],
            [self.intrusionDepth,"Intrusion Depth"]
        ]

    def __str__(self):
        out = ""
        out += "SMinPressure: " + str(self.SMinPressure) + "\n"
        out += "MLTFITPressure: " + str(self.MLTFITPressure) + "\n"
        out += "DThresholdPressure: " + str(self.DThresholdPressure) + "\n"
        out += "SGradientMaxPressure: " + str(self.SGradientMaxPressure) + "\n"
        out += "intrusionDepthPressure: " + str(self.intrusionDepthPressure) + "\n"
        out += "MLD Depth: " + str(self.foundMLD) + "\n"
        return out
