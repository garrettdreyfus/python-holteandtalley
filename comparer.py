import pickle
import math
import numpy as np
matlabOutput={}
pythonOutput={}

with open("pyOutput.pickle","rb") as f:
    pythonOutput = pickle.load(f)

with open("matOutput.pickle","rb") as f:
    matlabOutput = pickle.load(f)

deltas = {"salinityAlgoD":[],"densityAlgoD":[],"temperatureAlgoD":[],
        "densityThresholdD":[],"temperatureThresholdD":[],
        "densityGradientD":[],"temperatureGradientD":[],"tempMLTFITD":[],
        "debug":[],"id":[],"tempMLTFITIndex":[],"steepest":[],"densityMLTFITD":[],
        "salinityMLTFITD":[]
}

higherror = []
for i in pythonOutput:
    try:
        matout = matlabOutput[(int(i["platformNumber"]),int(i["cycleNumber"]))]
    except:
        print("could not find")
    ds = round(i["salinityAlgo"] - matout["salinityAlgo"],1)
    dt = round(i["tempAlgo"] - matout["tempAlgo"],1)
    dd = round(i["densityAlgo"] - matout["densityAlgo"],1)
    deltas["densityThresholdD"].append(
        round(i["densityThreshold"] - matout["densityThreshold"],1)
    )
    deltas["densityGradientD"].append(
        i["densityGradient"] - matout["densityGradient"]
    )
    deltas["temperatureThresholdD"].append(
        round(i["tempThreshold"] - matout["tempThreshold"],1)
    )
    deltas["temperatureGradientD"].append(
        i["tempGradient"] - matout["tempGradient"]
    )
    deltas["tempMLTFITD"].append(
        round(i["tempMLTFIT"] - matout["tempMLTFIT"],1)
    )
    #deltas["tempMLTFITIndex"].append(
        #i["tempMLTFITIndex"] - matout["tempMLTFITIndex"]
    #)
    deltas["steepest"].append(
        i["steepest"] - matout["steepest"]
    )
    deltas["densityMLTFITD"].append(
        round(i["densityMLTFIT"] - matout["densityMLTFIT"],1)
    )
    deltas["salinityMLTFITD"].append(
        round(i["salinityMLTFIT"] - matout["salinityMLTFIT"],1)
    )
    deltas["salinityAlgoD"].append(ds)
    deltas["temperatureAlgoD"].append(dt)
    deltas["densityAlgoD"].append(dd)
    deltas["debug"].append(i["debug"])
    #deltas["id"].append((int(i["platformNumber"]),i["cycleNumber"]))
    if abs(i["densityAlgo"] - matout["densityAlgo"]) > 10:
        higherror.append(((int(i["platformNumber"]),i["cycleNumber"]),i["tempAlgo"],matout["tempAlgo"],i["salinityAlgo"],matout["salinityAlgo"],i["debug"]))
#print(deltas)
#print(higherror)
print("="*20)
print("Salinities")
print("="*20)
print("mean sal d",sum(np.abs(deltas["salinityAlgoD"]))/len(deltas["salinityAlgoD"]))
print("max sal d",max(np.abs(deltas["salinityAlgoD"])))
print("sal stdev",(np.std(deltas["salinityAlgoD"])))
print("="*20)
print("Densities")
print("="*20)
print("mean density d",sum(np.abs(deltas["densityAlgoD"]))/len(deltas["densityAlgoD"]))
print("max density d",max(np.abs(deltas["densityAlgoD"])))
print("density stdev",(np.std(deltas["densityAlgoD"])))
print("="*20)
print("Temperatures")
print("="*20)
print("mean temp d",sum(np.abs(deltas["temperatureAlgoD"]))/len(deltas["temperatureAlgoD"]))
print("max temp d",max(np.abs(deltas["temperatureAlgoD"])))
print("temp stdev",(np.std(deltas["temperatureAlgoD"])))
print("="*20)
