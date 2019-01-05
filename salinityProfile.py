import numpy as np
class salinityProfile:
    def __init__(self,pressures=None,salinities=None,densities=None):
        self.mltfitline = []
        self.thermoclinefitline = []
        ##fnd reference pressure, this is done in holte and talley 
        ##supplementary matlab file
        startindex = np.argmin((np.asarray(pressures)-10)**2)
        self.salinities = salinities[startindex:]
        self.pressures = pressures[startindex:]
        self.densities = densities[startindex:]
        self.salinityGradients= self.generateGradientList()
        self.SMin = int(self.calculateSmin())
        self.SMinPressure = pressures[self.Smin]
        self.MLTFIT = int(self.calculateMLTFIT())
        self.MLTFITPressure = pressures[self.MLTFIT]
        self.DThreshold = int(self.calculateDThreshold())
        self.DThresholdPressure = pressures[self.DThreshold]
        #self.dT = self.calculateDeltaT()
        self.SGradientMax = int(self.calculateSGradientMax())
        self.SGradientMaxPressure = self.pressures[self.SGradientMax]
        self.intrusionDepth = int(self.calculateIntrusionDepth())
        self.intrusionPressure = self.pressures[self.intrusionDepth]
        self.foundMLD=0
        self.path=""
        #range from paper
        self.range = 25
        return

    ##find nearest pressure
    def findNearestPressureIndex(self,value):
        return (np.abs(np.asarray(self.pressures)- value)).argmin()

    def generateGradientList(self):
        tGS=[]
        for index in range(len(self.salinities)-1):
            dt = float(self.salinities[index] - self.salinities[index+1])
            dp = float(self.pressures[index] - self.pressures[index+1])
            tGS.append(dt/dp)
        smoothed=[0]*len(tGS)
        smoothed[0] = (tGS[0] + tGS[1])/2.0
        smoothed[-1] = (tGS[-1] + tGS[-2])/2.0
        for i in range(1,len(tGS)-1):
            smoothed[i] = (tGS[i-1]+tGS[i]+tGS[i+1])/3.0
        return tGS

    #The salinity minimum
    #Closest paper equivalent TMAX
    def calculateSMin(self):
        maxIndex = 0
        for i in range(len(self.salinities)):
            print(self.salinities[i])
            if self.salinities[i] >= self.salinities[maxIndex]:
                maxIndex = i
        return maxIndex

    #calculates the TTMLD or the salinity threshold mixed layer estimate
    # based on a salinity threshold of 0.03 from Boyer Montegut
    # Matlab equivalent: mldepthdens
    #closest paper equivalent: TTMLD
    def calculateDThreshold(self):
        for index in range(0,len(self.pressures)):
            if abs(self.salinities[index] - self.salinities[0]) > 0.03:
                return index
        return 0

    # calculates the MLTFIT or the intersection of the best fit mixed layer line
    # and the best fit of the thermocline
    #Matlab equivalent: upperdsmax
    #closest paper equivalent: MLTFIT
    def calculateMLTFIT(self):
        #Calculate the best fit line of the mixed layer
        errors = []
        fits = []
        #iterate through and polyfit over progressively increasing points
        for num in range(2,len(self.pressures)):
            out = np.polyfit(self.pressures[0:num],self.salinities[0:num],1,full=True)
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
        steepest = np.argmax(np.abs(self.salinityGradients))
        thermoclineFit = [self.salinityGradients[steepest],
            self.salinities[steepest]-self.salinityGradients[steepest]*self.pressures[steepest]
        ]
        self.thermoclinefitline = thermoclineFit
        depth = abs(float(thermoclineFit[1] - mltBestFit[1])/float(thermoclineFit[0] - mltBestFit[0]))
        return self.findNearestPressureIndex(depth)

    # The salinity difference across the mltfit or T(i mltfit) - T(i mltfit + 2 )
    def calculateDeltaT(self):
        if self.MLTFIT < len(self.salinities)-2:
            print(self.salinities[self.MLTFIT])
            print(self.salinities[self.MLTFIT+2])
            return float(self.salinities[self.MLTFIT] - self.salinities[self.MLTFIT+2])
        return len(self.salinities) -1

    #The salinity gradient threshold or max if threshold not met
    # Matlab: dsmin (confusingly named I know)
    def calculateSGradientMax(self):
        return np.argmax(self.salinityGradients)

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
            return self.Smin

    #Based on figure 9 from Holte and Talley 
    def mldWinterProfile(self):
        return

    #Based on figure 8 from Holte and Talley 
    def mldSummerProfile(self):
        return

    def findTemperatureMLD(self):
        if self.dT > 0.5 or self.dT < -0.25:
            self.foundMLD = self.mldSummerProfile()
        else:
            self.foundMLD =  self.mldWinterProfile()
        return self.foundMLD

    def __str__(self):
        out = ""
        out += "SMinPressure" + str() + "\n"
        out += "MLTFITPressure" + str() + "\n"
        out += "DThresholdPressure" + str() + "\n"
        out += "SGradientMaxPressure" + str() + "\n"
        out += "IntrusionDepthPressure" + str() + "\n"
        out += "" + str() + "\n"
        out += "" + str() + "\n"
        out += "" + str() + "\n"
        out += "" + str() + "\n"
        return out
