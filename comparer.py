import pickle
matlabOutput={}
pythonOutput={}

with open("pyOutput.pickle","rb") as f:
    pythonOutput = pickle.load(f)

with open("matOutput.pickle","rb") as f:
    matlabOutput = pickle.load(f)

deltas = {"salinityAlgoD":[],"densityAlgoD":[],"temperatureAlgoD":[],
        "densityThresholdD":[],"temperatureThresholdD":[],
        "densityGradientD":[],"temperatureGradientD":[],
        "debug":[]
}
for i in pythonOutput:
    try:
        matout = matlabOutput[(int(i["platformNumber"]),i["cycleNumber"])]
    except:
        print("could not find")
    ds = i["salinityAlgo"] - matout["salinityAlgo"]
    dt = i["tempAlgo"] - matout["tempAlgo"]
    dd = i["densityAlgo"] - matout["densityAlgo"]
    deltas["densityThresholdD"].append(
        i["densityThreshold"] - matout["densityThreshold"]
    )
    deltas["densityGradientD"].append(
        i["densityGradient"] - matout["densityGradient"]
    )
    deltas["temperatureThresholdD"].append(
        i["tempThreshold"] - matout["tempThreshold"]
    )
    deltas["temperatureGradientD"].append(
        i["tempGradient"] - matout["tempGradient"]
    )
    deltas["salinityAlgoD"].append(ds)
    deltas["temperatureAlgoD"].append(dt)
    deltas["densityAlgoD"].append(dd)
    deltas["debug"].append(i["debug"])
print(deltas)
