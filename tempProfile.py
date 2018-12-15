import numpy as np
class tempProfile:
    def __init__(self,temperatures=None,pressures=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        self.temperatures = temperatures
        self.pressures = pressures
        self.temperatureGradients= self.generateGradientList()
        self.TMax = int(self.calculateTMax())
        self.TTMLD = int(self.calculateTTMLD())
        self.MLTFIT = int(self.calculateMLTFIT())
        self.dT = int(self.calculateDeltaT())
        self.DTM = int(self.calculateDTM())
        self.TDTM = int(self.calculateTDTM())
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
        return tGS

    #The temperature maximum
    def calculateTMax(self):
        return np.argmax(self.temperatures)

    #calculates the TTMLD or the temperature threshold mixed layer estimate
    # based on a temperature threshold of 0.2
    def calculateTTMLD(self):
        for index in range(1,len(self.pressures)):
            if abs(self.pressures[index] - self.pressures[0]) > 0.2:
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
        steepest = np.argmax(self.temperatureGradients)
        thermoclineFit = [self.temperatureGradients[steepest],self.temperatures[steepest] -self.temperatureGradients[steepest]*steepest]
        depth = abs(float(thermoclineFit[1] - mltBestFit[1])/float(thermoclineFit[0] - mltBestFit[0]))
        self.thermoclinefitline = thermoclineFit
        return self.findNearestPressureIndex(depth)

        # The temperature difference across the mltfit or T(i mltfit) - T(i mltfit + 2 )
    def calculateDeltaT(self):
        if self.MLTFIT < len(self.temperatures)-2:
            return self.temperatures[self.MLTFIT] - self.temperatures[self.MLTFIT+2]
        return len(self.temperatures) -1

    #The density gradient threshold or max if threshold not met
    def calculateDTM(self):
        out = np.argwhere(np.array(self.temperatureGradients) > 0.005)
        if len(out) > 0:
            return out[0]
        return np.argmax(self.temperatureGradients)

    #The minimum of the depth of the temperature gradiet maximum and the temperature maximum
    def calculateTDTM(self):
        return min(self.DTM,self.TMax)

    # point j in figure 9
    def mldWinterPointJ(self,MLD):
        self.path += "J"
        if abs(MLD - self.TTMLD) > self.range and MLD == 0:
            MLD = self.TMax        
            if self.TMax == 1:
                MLD = self.TTMLD
            if self.Tmax > self.TTMLD:
                MLD= self.TTMLD
            return MLD
        else:
            return MLD

    # point F in figure 9
    def mldWinterPointF(self,MLD):
        self.path += "F"
        if ((abs(self.TDM - self.MLTFIT) < self.range and abs(self.TDM-self.TTMLD)<self.range) or
            (abs(self.TDM - self.TTMLD) < self.range and abs(self.TTMLD - self.MLTFIT) < self.range) or
            (abs(self.TDM - self.MLTFIT) < self.range and abs(self.TTMLD - self.MLTFIT) < self.range)
           ):
            MLD = self.MLTFIT
        if MLD > self.TTMLD:
            MLD=self.TTMLD
        self.mldWinterPointJ(MLD)
        
    # point H in figure 9
    def mldWinterPointH(self,MLD):
        self.path += "H"
        if self.MLTFIT - self.TTMLD < self.range:
            MLD = self.MLTFIT
            return self.mldWinterPointJ(MLD)
        else:
            MLD = self.DTM
            if MLD>self.TTMLD:
                MLD = self.TTMLD
            return self.mldWinterPointJ(MLD)

    #Based on figure 9 from Holte and Talley 
    def mldWinterProfile(self):
        MLD = -99999
        if (abs(self.pressures[self.MLTFIT] - self.pressures[self.TTMLD]) < self.range and
                abs(self.TDTM - self.TTMLD) > self.range and
                self.MLTFIT <self.TDTM):
            MLD = self.MLTFIT
            return self.mldWinterPointJ(MLD)
        else:
            if self.TDTM > self.range:
                MLD = self.TDTM
                return self.mldWinterPointF(MLD)
            else:
                return self.mldWinterPointH(MLD)
        return

    #Based on figure 8 from Holte and Talley 
    def mldSummerProfile(self):
        MLD=self.MLTFIT
        if self.dT<0 and MLD > self.TTMLD:
            MLD = self.TTMLD
        if MLD > self.TTMLD:
            if self.TMax < self.TTMLD and self.TMax > self.range:
                MLD = self.TMAx
            else:
                MLD = self.TTMLD
        self.foundMLD = MLD
        return self.foundMLD

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
        return out
