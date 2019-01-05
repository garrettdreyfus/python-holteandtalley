from tempProfile import tempProfile
from densityProfile import densityProfile
from salinityProfile import salinityProfile
import numpy as np
import gsw as gsw
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import os

def profilePlotter(x,y,line,depths,variable,lat,lon,path):
    plt.plot(x,y)
    plt.ylim(500,0)
    plt.axhline(line)
    for i in depths:
        #Make random markers with labels
        print(i)
        plt.plot(x[i[0]],y[i[0]],'ro')
    plt.xlabel("Densities (kg/m^3)")
    plt.ylabel("Pressures (dbar)")
    plt.title(str("lat : "+lat+
                "lon : "+lon+
                "path :  "+path + "  mlt: " + path))
    plt.show()

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
        temps = tempsOut
        densities = densitiesOut
        pressuress = pressuresOut
        salinities = salts
        #b = densityProfile(densitiesOut,ps)
        s = salinityProfile(pressures,temps,salinities,densities)
        s.findMLD()
        d = densityProfile(pressures,temps,salts,densities,sp=s)
        d.findMLD()
        #print(s)
        #profilePlotter(salinities,pressures,s.foundMLD,s.importantDepths(),"Salinities",
            #str(dataset.variables["LATITUDE"][0]),
            #str(dataset.variables["LONGITUDE"][0]),
            #str(s.path)
            #)
        print(d)
        profilePlotter(densities,pressures,d.foundMLD,d.importantDepths(),"Densities",
            str(dataset.variables["LATITUDE"][0]),
            str(dataset.variables["LONGITUDE"][0]),
            str(d.path)
            )
