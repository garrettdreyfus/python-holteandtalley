from tempProfile import tempProfile
from densityProfile import densityProfile
import numpy as np
import gsw as gsw
from netCDF4 import Dataset
import matplotlib.pyplot as plt


import os
for file in os.listdir("profiles"):
    if file.endswith(".nc"):
        plt.figure()
        filename = os.path.join("profiles", file)
        dataset = Dataset(filename)
        pressures = dataset.variables["PRES"][:][0]
        salts = dataset.variables["PSAL"][:][0]
        temps = dataset.variables["TEMP"][:][0]
        tempsOut=[]
        pressuresOut=[]
        densitiesOut=[]
        for index in range(len(pressures)):
            if pressures[index] != "_":
                pres = pressures[index]
                psal = salts[index]
                temp = temps[index]
                temp = gsw.conversions.CT_from_t(psal,temp,pres)
                tempsOut.append(temp)
                pressuresOut.append(float(pres))
                densitiesOut.append(float(gsw.sigma0(psal,temp)))
        ts = tempsOut
        ps = pressuresOut
        #try:
        #except:
            #print(filename)
            #break
        #try:
        a = tempProfile(ts,ps)
        b = densityProfile(densitiesOut,ps)
        #except Exception as e:
            #print(e)
            #print("did not converge")
            #break
        plt.plot(ts,ps)
        #plt.ylim(max(ps),min(ps))
        plt.ylim(500,0)
        plt.plot(ts[a.TMax],ps[a.TMax],'bo')
        plt.plot(ts[a.TTMLD],ps[a.TTMLD],'gs')
        plt.plot(ts[a.MLTFIT],ps[a.MLTFIT],'o',color="orange")
        plt.plot(ts[a.DTM],ps[a.DTM],'g>')
        plt.plot(ts[a.TDTM],ps[a.TDTM],'bs')
        a.findTemperatureMLD()
        #plt.ylim(500,0)
        plt.axhline(a.foundMLD)
        print(a.temperatureGradients)
        plt.xlabel("Temperature (C)")
        plt.ylabel("Pressures (dbar)")
        plt.title(str("lat : "+str(dataset.variables["LATITUDE"][0])+
                    "lon : "+str(dataset.variables["LONGITUDE"][0])+
                    "path :  "+a.path + "  mlt: " + str(a.foundMLD)))
        plt.figure()
        plt.plot(densitiesOut,ps)
        #plt.ylim(max(ps),min(ps))
        plt.ylim(500,0)
        plt.plot(densitiesOut[b.TMax],ps[b.TMax],'bo')
        plt.plot(densitiesOut[b.TTMLD],ps[b.TTMLD],'gs')
        plt.plot(densitiesOut[b.MLTFIT],ps[b.MLTFIT],'o',color="orange")
        plt.plot(densitiesOut[b.DTM],ps[b.DTM],'g>')
        plt.plot(densitiesOut[b.TDTM],ps[b.TDTM],'bs')
        b.findTemperatureMLD()
        #plt.ylim(500,0)
        plt.axhline(b.foundMLD)
        print(b.densityGradients)
                
        plt.xlabel("Densities (kg/m^3)")
        plt.ylabel("Pressures (dbar)")
        plt.title(str("lat : "+str(dataset.variables["LATITUDE"][0])+
                    "lon : "+str(dataset.variables["LONGITUDE"][0])+
                    "path :  "+b.path + "  mlt: " + str(b.foundMLD)))
        print(a)
plt.show()
