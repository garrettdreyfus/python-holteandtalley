from tempProfile import tempProfile
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
        pressures = dataset.variables["pres"][:][0]
        salts = dataset.variables["psal"][:][0]
        temps = dataset.variables["temp"][:][0]
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
        try:
            a = tempProfile(ts,ps)
        except:
            print(filename)
        plt.plot(ts,ps)
        plt.ylim(max(ps),min(ps))
        plt.plot(ts[a.TMax],ps[a.TMax],'bo')
        plt.plot(ts[a.TTMLD],ps[a.TTMLD],'ro')
        plt.plot(ts[a.MLTFIT],ps[a.MLTFIT],'go')
        plt.plot(ts[a.dT],ps[a.dT],'bx')
        plt.plot(ts[a.DTM],ps[a.DTM],'gx')
        plt.plot(ts[a.TDTM],ps[a.TDTM],'rx')
        a.findTemperatureMLD()
        y_vals = []
        #for i in ps:
            #y_vals.append(i*a.mltfitline[0]+a.mltfitline[1])
        #plt.plot(y_vals,ps)
        #y_vals = []
        #for i in ps:
            #y_vals.append(i*a.thermoclinefitline[0]+a.thermoclinefitline[1])
        #plt.plot(y_vals,ps)
        #plt.ylim(500,0)
        plt.axhline(ps[a.foundMLD])
        plt.xlabel("Temperature (C)")
        plt.ylabel("Pressures (dbar)")
        plt.title(str("lat : "+str(dataset.variables["latitude"][0])+
                    "lon : "+str(dataset.variables["longitude"][0])+
                    "path :  "+a.path + "mlt index: " + str(a.foundMLD)))
        #print(a)
plt.show()
