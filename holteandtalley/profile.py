import numpy as np
class Profile():
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
            
    def generateGradientListNoSmooth(self,values):
        tGS=[]
        for index in range(len(values)-1):
            dt = float(values[index] - values[index+1])
            dp = float(self.pressures[index] - self.pressures[index+1])
            tGS.append(dt/dp)
        return tGS[1:-1]

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
                self.mltfitline=mltBestFit
                self.mltfitindex = index-1
                break
        #find thermoclineFit
        steepest = np.argwhere(np.abs(gradients) == np.max(np.abs(gradients)))[-1][0] +1

        thermoclineFit = np.polyfit(self.pressures[steepest-1:steepest+2],values[steepest-1:steepest+2],1,full=True)[0]
        if thermoclineFit[0] != mltBestFit[0]:
            depth = abs(float(thermoclineFit[1] - mltBestFit[1])/float(thermoclineFit[0] - mltBestFit[0]))
        else:
            depth = self.pressures[steepest]

        ## a switch to choose whether to use nearest pressure or not
        if False:
            return self.findNearestPressureIndex(depth), depth
        else:
            return self.findNearestPressureIndex(depth),self.pressures[self.findNearestPressureIndex(depth)]# depth

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
            return 0 
        else:
            return 1
    #The salinity gradient threshold or max if threshold not met
    # Matlab: dsmin (confusingly named I know)
    def calculateSGradientMax(self):
        return np.argwhere(np.abs(self.salinityGradients) == np.max(np.abs(self.salinityGradients)))[-1][0]+1





