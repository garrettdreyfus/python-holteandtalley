from .profile import Profile
import numpy as np
class tempProfile(Profile):
    def __init__(self,pressures,temperatures):
        self.mltfitline = []
        self.thermoclinefitline = []
        self.debug = 0
        ##fnd reference pressure, this is done in holte and talley 
        ##supplementary matlab file
        startindex = np.argmin((np.asarray(pressures)-10.0)**2)
        self.temperatures = temperatures[startindex:]
        self.pressures = pressures[startindex:]
        self.temperatureGradients= self.generateGradientListNoSmooth(self.temperatures)
        self.TMax = int(self.calculateTMax())
        self.TMaxPressure = self.pressures[self.TMax]
        self.MLTFIT, self.MLTFITPressure = self.calculateMLTFIT(self.temperatures,self.temperatureGradients)
        self.TTMLD = int(self.calculateTTMLD())
        self.TTMLDPressure = self.interpolateTTMLD()
        self.dT = self.calculateDeltaT()
        self.DTM = int(self.calculateDTM())
        self.DTMPressure = self.pressures[self.DTM]
        self.TDTM = int(self.calculateTDTM())
        self.TDTMPressure = self.pressures[self.TDTM]
        self.foundMLD=0
        self.path=""
        #range from paper
        self.range = 25
        return

    #The temperature maximum
    def calculateTMax(self):
        return np.where(self.temperatures == np.max(self.temperatures))[-1][-1]

    #calculates the TTMLD or the temperature threshold mixed layer estimate
    # based on a temperature threshold of 0.2
    def calculateTTMLD(self):
        for index in range(0,len(self.pressures)):
            if abs(self.temperatures[index] - self.temperatures[0]) > 0.2:
                return index
        return 0

    def interpolateTTMLD(self):
        if self.TTMLD:
            thresholdIndex = self.TTMLD
        else:
            thresholdIndex = int(self.calculateTTMLD())
        if thresholdIndex ==0:
            return self.pressures[thresholdIndex]
        deltaP = (self.pressures[thresholdIndex] - self.pressures[thresholdIndex-1])
        deltaT = (self.temperatures[thresholdIndex] - self.temperatures[thresholdIndex-1])
        if self.temperatures[thresholdIndex] < self.temperatures[thresholdIndex-1]:
            return self.pressures[thresholdIndex-1] + (deltaP/deltaT)*((self.temperatures[0]-0.2)-self.temperatures[thresholdIndex-1])
        elif self.temperatures[thresholdIndex] > self.temperatures[thresholdIndex-1]:
            return self.pressures[thresholdIndex-1] + (deltaP/deltaT)*((self.temperatures[0]+0.2)-self.temperatures[thresholdIndex-1])
        #return self.pressures[thresholdIndex]# + (deltaP/deltaT)*(0.03-abs(self.temperatures[thresholdIndex-1]))

    # The temperature difference across the mltfit or T(i mltfit) - T(i mltfit + 2 )
    def calculateDeltaT(self):
        if self.MLTFIT < len(self.temperatures)-2:
            return float(self.temperatures[self.MLTFIT] - self.temperatures[self.MLTFIT+2])
        return len(self.temperatures) -1

    #The density gradient threshold or max if threshold not met
    def calculateDTM(self):
        maxIndex=0
        for i in range(len(self.temperatureGradients)):
            if abs(self.temperatureGradients[i]) > 0.005 :
                return i+1
            elif abs(self.temperatureGradients[i]) > self.temperatureGradients[maxIndex]:
                maxIndex=i+1
        return maxIndex
    #The minimum of the depth of the temperature gradiet maximum and the temperature maximum
    def calculateTDTM(self):
        steepest = np.where(self.temperatureGradients == np.max(self.temperatureGradients))[-1][-1] + 1
        if self.pressures[steepest] < self.TMaxPressure:
            return steepest
        else:
            return self.TMax

    # point j in figure 9
    def mldWinterPointJ(self,MLD):
        self.path += "J"
        if abs(MLD - self.TTMLDPressure) > self.range and MLD == 0:
            MLD = self.TMaxPressure        
            self.debug="TMAX"
            if self.TMaxPressure == 1:
                MLD = self.TTMLDPressure
                self.debug="TTMLD"
            if self.TMaxPressure > self.TTMLDPressure:
                self.debug="TTMLD"
                MLD= self.TTMLDPressure
            return MLD
        else:
            return MLD

    # point F in figure 9
    def mldWinterPointF(self,MLD):
        self.path += "F"
        if ((abs(self.TDTMPressure - self.TDTMPressure) < self.range and abs(self.TDTMPressure-self.TTMLDPressure)<self.range) or
            (abs(self.TDTMPressure - self.TDTMPressure) < self.range and abs(self.TTMLDPressure - self.MLTFITPressure) < self.range) or
            (abs(self.TDTMPressure - self.TDTMPressure) < self.range and abs(self.TTMLDPressure - self.MLTFITPressure) < self.range)
           ):
            MLD = self.MLTFITPressure
            self.debug="MLTFIT zap"
        if MLD > self.TTMLDPressure:
            MLD=self.TTMLDPressure
            self.debug="TTMLD"
        return self.mldWinterPointJ(MLD)
        
    # point H in figure 9
    def mldWinterPointH(self,MLD):
        self.path += "H"
        if (self.MLTFITPressure - self.TTMLDPressure) < self.range:
            MLD = self.MLTFITPressure
            self.debug="MLTFIT zip"
            return self.mldWinterPointJ(MLD)
        else:
            MLD = self.DTMPressure
            self.debug="DTM"
            if MLD>self.TTMLDPressure:
                MLD = self.TTMLDPressure
                self.debug="TTMLD"
            return self.mldWinterPointJ(MLD)

    #Based on figure 9 from Holte and Talley 
    def mldWinterProfile(self):
        MLD = -99999
        if (abs(self.MLTFITPressure - self.TTMLDPressure) < self.range and
                abs(self.TDTMPressure - self.TTMLDPressure) > self.range and
                self.MLTFITPressure <self.TDTMPressure):
            MLD = self.MLTFITPressure
            self.debug="MLTFIT zop"
            return self.mldWinterPointJ(MLD)
        else:
            if self.TDTMPressure > self.pressures[0] + self.range:
                MLD = self.TDTMPressure
                self.debug="TDTM"
                return self.mldWinterPointF(MLD)
            else:
                return self.mldWinterPointH(MLD)

    #Based on figure 8 from Holte and Talley 
    def mldSummerProfile(self):
        MLD=self.MLTFITPressure
        self.debug="MLTFIT no action"
        if self.dT<0 and MLD > self.TTMLDPressure:
            MLD = self.TTMLDPressure
            self.debug="TTMLD"
        if MLD > self.TTMLDPressure:
            if self.TMaxPressure < self.TTMLDPressure and self.TMaxPressure > self.range:
                MLD = self.TMaxPressure
                self.debug="TMax"
            else:
                MLD = self.TTMLDPressure
                self.debug="TTMLD"
        return MLD

    def findMLD(self):
        if self.dT > 0.5 or self.dT < -0.25:
            self.season = 0
            self.foundMLD = self.mldSummerProfile()
        else:
            self.season = 1
            self.foundMLD =  self.mldWinterProfile()
        return self.foundMLD
    def __str__(self):
        out = ""
        out += "Tmax: " + str(self.TMax) + "\n"
        out += "TTMLD: " + str(self.TTMLD) + "\n"
        out += "MLTFIT: " + str(self.MLTFIT) + "\n"
        out += "MLTFIT Line: " + str(self.mltfitline) + "\n"
        out += "Thermocline Line: " + str(self.thermoclinefitline) + "\n"
        out += "dT: " + str(self.dT) + "\n"
        out += "DTM: " + str(self.DTM) + "\n"
        out += "TDDTM: " + str(self.TDTM) + "\n"
        out += "Found MLD: " + str(self.foundMLD) + "\n"
        return out
