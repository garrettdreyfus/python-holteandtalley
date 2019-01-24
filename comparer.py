import pickle
matlabOutput={}
pythonOutput={}

with open("pyOutput.pickle","rb") as f:
    pythonOutput = pickle.load(f)

with open("matOutput.pickle","rb") as f:
    matlabOutput = pickle.load(f)

deltas = {"salinityAlgoD":[],"densityAlgoD":[],"temperatureAlgoD":[],
        "densityThresholdD":[],"temperatureThresholdD":[],
        "densityGradientD":[],"temperatureGradientD":[],"tempMLTFITD":[],
        "debug":[],"id":[],"mltFitIndex":[],"steepest":[]
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
    deltas["mltFitIndex"].append(
        i["mltFitIndex"] - matout["mltFitIndex"]
    )
    deltas["steepest"].append(
        i["steepest"] - matout["steepest"]
    )
    deltas["salinityAlgoD"].append(ds)
    deltas["temperatureAlgoD"].append(dt)
    deltas["densityAlgoD"].append(dd)
    deltas["debug"].append(i["debug"])
    #deltas["id"].append((int(i["platformNumber"]),i["cycleNumber"]))
    if abs(dt) > 10:
        higherror.append(((int(i["platformNumber"]),i["cycleNumber"]),i["tempAlgo"],matout["tempAlgo"],i["debug"]))

print(deltas)
print(higherror)
