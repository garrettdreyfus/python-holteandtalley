import numpy as np
class tempProfile:
    def __init__(self,temperatures=None,pressures=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        self.temperatures = temperatures
        self.pressures = pressures
        self.temperatureGradients= self.generateGradientList()
        self.TMax = int(self.calculateTMax())
        self.TMaxPressure = pressures[self.TMax]
        self.MLTFIT = int(self.calculateMLTFIT())
        self.MLTFITPressure = pressures[self.MLTFIT]
        self.TTMLD = int(self.calculateTTMLD())
        self.TTMLDPressure = pressures[self.TTMLD]
        self.dT = self.calculateDeltaT()
        self.DTM = int(self.calculateDTM())
        self.DTMPressure = self.pressures[self.DTM]
        self.TDTM = int(self.calculateTDTM())
        self.TDTMPressure = self.pressures[self.TDTM]
        self.foundMLD=0
        self.path=""
        self.path=""
        #range from paper
        self.range = 25
        return

    ##find nearest pressure
    def findNearestPressureIndex(self,value):
        return (np.abs(np.asarray(self.pressures)- value)).argmin()

    def generateGradientList(self):
        tGS=[]
        for index in range(len(self.temperatures)-1):
            dt = float(self.temperatures[index] - self.temperatures[index+1])
            dp = float(self.pressures[index] - self.pressures[index+1])
            tGS.append(dt/dp)
        smoothed=[0]*len(tGS)
        smoothed[0] = (tGS[0] + tGS[1])/2.0
        smoothed[-1] = (tGS[-1] + tGS[-2])/2.0
        for i in range(1,len(tGS)-1):
            smoothed[i] = (tGS[i-1]+tGS[i]+tGS[i+1])/3.0
        return tGS

    #The temperature maximum
    def calculateTMax(self):
        return np.argmax(self.temperatures)

    #calculates the TTMLD or the temperature threshold mixed layer estimate
    # based on a temperature threshold of 0.2
    def calculateTTMLD(self):
        ref = np.argmin(abs(np.asarray(self.pressures) -10))
        print(ref)
        for index in range(0,len(self.pressures)):
            if abs(self.temperatures[index] - self.temperatures[ref]) > 0.2:
                return index
        return 0

    # calculates the MLTFIT or the intersection of the best fit mixed layer line
    # and the best fit of the thermocline
    def calculateMLTFIT(self):
        #Calculate the best fit line of the mixed layer
        errors = []
        fits = []
        #iterate through and polyfit over progressively increasing points
        for num in range(2,len(self.pressures)):
            out = np.polyfit(self.pressures[0:num],self.temperatures[0:num],1,full=True)
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
        steepest = np.argmax(np.abs(self.temperatureGradients))
        thermoclineFit = [self.temperatureGradients[steepest],
            self.temperatures[steepest]-self.temperatureGradients[steepest]*self.pressures[steepest]
        ]
        self.thermoclinefitline = thermoclineFit
        depth = abs(float(thermoclineFit[1] - mltBestFit[1])/float(thermoclineFit[0] - mltBestFit[0]))
        return self.findNearestPressureIndex(depth)

        # The temperature difference across the mltfit or T(i mltfit) - T(i mltfit + 2 )
    def calculateDeltaT(self):
        if self.MLTFIT < len(self.temperatures)-2:
            print(self.temperatures[self.MLTFIT])
            print(self.temperatures[self.MLTFIT+2])
            return float(self.temperatures[self.MLTFIT] - self.temperatures[self.MLTFIT+2])
        return len(self.temperatures) -1

    #The density gradient threshold or max if threshold not met
    def calculateDTM(self):
        maxIndex=0
        for i in range(len(self.temperatureGradients)):
            if abs(self.temperatureGradients[i]) > 0.005:
                return i
            elif abs(self.temperatureGradients[i]) > self.temperatureGradients[maxIndex]:
                maxIndex=i
        return maxIndex

    #The minimum of the depth of the temperature gradiet maximum and the temperature maximum
    def calculateTDTM(self):
        steepest=0
        for i in range(len(self.temperatureGradients)):
            if  (self.temperatureGradients[i]) > self.temperatureGradients[steepest]:
                steepest=i
        if self.pressures[steepest] < self.pressures[self.TMax]:
            return steepest
        else:
            return self.TMax

    # point j in figure 9
    def mldWinterPointJ(self,MLD):
        self.path += "J"
        if abs(MLD - self.TTMLDPressure) > self.range and MLD == 0:
            MLD = self.TMaxPressure        
            if self.TMaxPressure == 1:
                MLD = self.TTMLDPressure
            if self.TmaxPressure > self.TTMLDPressure:
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
        if MLD > self.TTMLDPressure:
            MLD=self.TTMLDPressure
        return self.mldWinterPointJ(MLD)
        
    # point H in figure 9
    def mldWinterPointH(self,MLD):
        self.path += "H"
        if self.MLTFITPressure - self.TTMLDPressure < self.range:
            MLD = self.MLTFITPressure
            return self.mldWinterPointJ(MLD)
        else:
            MLD = self.DTMPressure
            if MLD>self.TTMLDPressure:
                MLD = self.TTMLDPressure
            return self.mldWinterPointJ(MLD)

    #Based on figure 9 from Holte and Talley 
    def mldWinterProfile(self):
        MLD = -99999
        if (abs(self.MLTFITPressure - self.TTMLDPressure) < self.range and
                abs(self.TDTMPressure - self.TTMLDPressure) > self.range and
                self.MLTFITPressure <self.TDTMPressure):
            MLD = self.MLTFITPressure
            return self.mldWinterPointJ(MLD)
        else:
            if self.TDTMPressure > self.range:
                MLD = self.TDTMPressure
                return self.mldWinterPointF(MLD)
            else:
                return self.mldWinterPointH(MLD)

    #Based on figure 8 from Holte and Talley 
    def mldSummerProfile(self):
        MLD=self.MLTFITPressure
        if self.dT<0 and MLD > self.TTMLDPressure:
            MLD = self.TTMLDPressure
        if MLD > self.TTMLDPressure:
            if self.TMaxPressure < self.TTMLDPressure and self.TMaxPressure > self.range:
                MLD = self.TMaxPressure
            else:
                MLD = self.TTMLDPressure
        return MLD

    def findTemperatureMLD(self):
        if self.dT > 0.5 or self.dT < -0.25:
            self.foundMLD = self.mldSummerProfile()
        else:
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
