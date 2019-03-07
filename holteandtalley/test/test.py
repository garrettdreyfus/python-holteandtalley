from holteandtalley import HolteAndTalley
import numpy as np
import gsw as gsw
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import os
import pickle

def profilePlotter(x,y,line,depths,variable,lat,lon,path):
    plt.plot(x,y)
    plt.ylim(500,0)
    plt.axhline(line)
    for i in depths:
        #Make random markers with labels
        plt.plot(x[i[0]],y[i[0]],'ro')
    plt.xlabel("Densities (kg/m^3)")
    plt.ylabel("Pressures (dbar)")
    plt.title(str("lat : "+lat+
                "lon : "+lon+
                "path :  "+path + "  mlt: " + path))
    plt.show()
out = []
for file in os.listdir("profiles"):
    if file.endswith(".nc"):
        filename = os.path.join("profiles", file)
        dataset = Dataset(filename)
        cycleNumber = dataset.variables["CYCLE_NUMBER"][:][0]
        pressures = dataset.variables["PRES"][:][0]
        salts = dataset.variables["PSAL"][:][0]
        temps = dataset.variables["TEMP"][:][0]

        pressuresa = dataset.variables["PRES_ADJUSTED"][:][0]
        saltsa = dataset.variables["PSAL_ADJUSTED"][:][0]
        tempsa = dataset.variables["TEMP_ADJUSTED"][:][0]

        if pressuresa[0] <99999:
            pressures =pressuresa
            salts = saltsa
            temps = tempsa
        if len(np.where(np.where(pressures) ==99999))>1 and np.where(np.where(pressures) ==99999)[0] <5:
            print("well shucks")
        else:
            tempsOut=[]
            pressuresOut=[]
            densitiesOut=[]
            salinitiesOut=[]
            for index in range(len(pressures)):
                if pressures[index] != "_":
                    pres = pressures[index]
                    psal = salts[index]
                    temp = temps[index]
                    temp = gsw.conversions.CT_from_t(psal,temp,pres)
                    tempsOut.append(temp)
                    pressuresOut.append(float(pres))
                    densitiesOut.append(float(gsw.sigma0(psal,temp)))
                    salinitiesOut.append(psal)
            #temps = tempsOut
            densities = densitiesOut
            pressuress = pressuresOut
            salinities = salts
            #b = densityProfile(densitiesOut,ps)
            line = {}
            h = HolteAndTalley(pressures,temps,salinities,densities)
            t = h.temp
            line["tempAlgo"] = t.findMLD()
            num = ""
            for i in dataset.variables["PLATFORM_NUMBER"][:][0]:
                if len(str(i)) >= 4:
                    num+=str(i)[2]
            line["platformNumber"] = num 
            s = h.salinity
            line["salinityAlgo"] = s.findMLD()
            d=h.density
            line["densityAlgo"] = d.findMLD()
            line["tempThreshold"] = t.TTMLDPressure
            line["tempGradient"] = t.DTMPressure
            line["densityThreshold"] = d.DThresholdPressure 
            line["densityGradient"] = d.DGradientThresholdPressure
            line["cycleNumber"] = cycleNumber
            line["tempMLTFIT"] = t.MLTFITPressure
            line["tempMLTFITIndex"] = t.mltfitindex
            line["densityMLTFIT"] = d.MLTFITDensityPressure
            line["salinityMLTFIT"] = s.MLTFITSalinityPressure
            line["salinityMLTFITIndex"] = s.mltfitindex
            line["debug"] = d.debug
            line["steepest"] = s.SGradientMax
            out.append(line)
            ##if int(num) == 3900623 and cycleNumber == 3:
                #print("MLTFIT: ", s.MLTFITPressure)
                ##print("MLTFITINDEX: ",s.mltfitindex)
                #print("D Threshold: ",s.DThresholdPressure)
                ##print("DMIN: ",s.SGradientMax)
                ##print(s.intrusionDepthPressure)
                ##print(d.MLDT,d.DThresholdPressure)
                ##=====
                ##print("MLTFIT: ", d.MLTFITPressure)
                ##print("D Threshold: ",d.DThresholdPressure)
                ##print("D Gradient Threshold: ",d.DGradientThresholdPressure)
                ##print("DMIN: ",d.DMinPressure)
                ##====
                ##print("TMAX: ",d.TMaxPressure)
                ##print(d.debug)
                ##print(t.mltfitindex)
                ##print(t.temperatureGradients)
                #profilePlotter(s.salinities,s.pressures,s.MLTFITPressure,[],"","","","")
with open("pyOutput.pickle","wb") as f:
    pickle.dump(out,f)
