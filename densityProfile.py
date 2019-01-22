import numpy as np
from tempProfile import tempProfile
from salinityProfile import tempProfile
class densityProfile:
    def __init__(self,pressures,temperatures,salinities,densities,tp=None,sp=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        if tp == None:
            self.tp = tempProfile(pressures,temperatures)
        else:
            self.tp = tp
        if sp == None:
            self.sp = tempProfile(pressures,temperatures)
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
        self.MLTFIT = int(self.calculateMLTFIT(self.densities,self.densityGradients))
        self.MLTFITPressure = pressures[self.MLTFIT]
        self.DThreshold = int(self.calculateDThreshold())
        self.DThresholdPressure = pressures[self.DThreshold]
        self.DGradientThreshold = int(self.calculateDGradientThreshold())
        self.DGradientThresholdPressure = pressures[self.DGradientThreshold]
        self.densityTest = sp.calculateDensityTest()
        #self.SGradientMax = int(self.calculateSGradientMax())
        #self.SGradientMaxPressure = self.pressures[self.SGradientMax]
        #self.intrusionDepth = int(self.calculateIntrusionDepth())
        #self.intrusionDepthPressure = self.pressures[self.intrusionDepth]
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
        smoothed=[0]*len(tGS)
        smoothed[0] = (tGS[0] + tGS[1])/2.0
        smoothed[-1] = (tGS[-1] + tGS[-2])/2.0
        for i in range(1,len(tGS)-1):
            smoothed[i] = (tGS[i-1]+tGS[i]+tGS[i+1])/3.0
        return smoothed

    #The Density Maximum
    #Closest paper equivalent TMAX
    def calculateDMin(self):
        return np.argmin(self.densities)

    #calculates the TTMLD or the salinity threshold mixed layer estimate
    # based on a salinity threshold of 0.03 from Boyer Montegut
    # Matlab equivalent: mldepthdens
    #closest paper equivalent: TTMLD
    def calculateDThreshold(self):
        for index in range(0,len(self.pressures)):
            if abs(self.densities[index] - self.densities[0]) > 0.03:
                return index
        return 0

    #Calculates the densityGradient threshold or max if criteria not met
    #closest to DTM
    def calculateDGradientThreshold(self):
        for index in range(0,len(self.densityGradients)):
            if abs(self.densityGradients[index]) > 0.0005:
                return index+1
        return np.argmax(self.densityGradients)+1

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
            ddiff = self.densities(densityGradientMax-1) - self.densities(densityGradientMax+1)
        #various constants from paper
        if ddiff > -0.06 and self.tp.dT > 0.5:
            return True
        if ddiff > -0.06 and self.tp.dT < -0.25:
            return False
        if self.tp.dT > -0.25 and self.tp.dT < 0.5:
            return True 
        else:
            return False

    #The salinity gradient threshold or max if threshold not met
    # Matlab: dsmin (confusingly named I know)
    def calculateSGradientMax(self):
        return np.argmax(self.salinityGradients)

    #The minimum of the depth of the salinity gradiet maximum and the salinity minimum
    #In matlab file: dtandtmin
    def calculateIntrusionDepth(self):
        steepest=0
        for i in range(len(self.densityGradients)):
            if  (self.densityGradients[i]) > self.densityGradients[steepest]:
                steepest=i
        if self.pressures[steepest] < self.pressures[self.DMax]:
            return steepest
        else:
            return self.DMax

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldWinterProfile(self):
        MLD = self.DThresholdPressure
        if self.tp.TTMLDPressure < MLD:
            MLD = self.tp.TTMLDPressure
        if self.MLTFITPressure < self.DThresholdPressure and self.MLTFITPressure > self.range:
            MLD = self.MLTFITPressure
        if self.tp.TDTMPressure > self.range and self.tp.TDTMPressure < self.DThresholdPressure:
            MLD= self.tp.TDTMPressure
            if abs(self.tp.TMaxPressure - self.DThresholdPressure) < abs(self.tp.TDTMPressure - self.MLTFITPressure):
                MLD = self.tp.TMaxPressure
            if abs(self.MLDS - self.DThresholdPressure) < self.range and self.MLDS < self.DThresholdPressure:
                MLD = min(self.DThresholdPressure,self.MLDS)
        if abs(self.MLDT - self.MLDS) < self.range and abs(min(self.MLDT,self.MLDS))-MLD > self.range:
            MLD = min(self.MLDT,self.MLDS)
        if MLD > self.DGradientThresholdPressure and abs(self.DGradientThresholdPressure - self.MLDT) < abs(MLD-self.MLDT):
            MLD = self.DGradientThresholdPressure
        if self.MLTFITPressure == self.sp.MLTFITPressure and abs(self.sp.MLTFITPressure - self.DThresholdPressure) < self.range:
            MLD= self.MLTFITPressure
        if self.MLDT == self.DMinPressure:
            MLD = self.DMinPressure
        return MLD

    #Based on find_mld.m from Holte and Talley Suplementary materials
    def mldSummerProfile(self):
        MLD = self.MLTFITPressure
        if MLD > self.MLTFITPressure:
            MLD = self.DThresholdPressure
        ## If any two of these three conditions are true
        conditions = [abs(self.MLDS-self.MLDT),abs(self.MLTFITPressure-self.MLDT),abs(self.MLDS - self.MLTFITPressure)]
        s =0
        for a in conditions:
            if a <self.range:
                s+=1
        if s >1:
            MLD = self.MLTFITPressure
        if abs(self.MLDS - self.DThresholdPressure) < self.range and self.MLDS != self.DThresholdPressure:
            if self.DThresholdPressure < self.MLDS:
                MLD = self.DThresholdPressure
            else:
                MLD = self.MLDS
            if self.MLTFITPressure == self.DThresholdPressure:
                MLD = self.MLTFITPressure
        if MLD > self.DGradientThresholdPressure and abs(self.DGradientThresholdPressure - self.MLDT):
           MLD = self.DGradientThresholdPressure 
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
