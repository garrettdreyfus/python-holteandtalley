from scipy.io import loadmat
import pickle
mldinfo =loadmat('mldinfo.mat')["mldinfo"]
out={}
print(mldinfo)
for i in mldinfo:
    line={}
    line["floatNumber"] = i[0]
    line["cycleNumber"] = i[26]
    line["tempMLTFIT"] = i[27]
    line["tempMLTFITIndex"] = i[28]
    line["densityMLTFIT"] = i[30]
    line["salinityMLTFIT"] = i[31]
    line["steepest"] = i[29]
    line["tempAlgo"] = i[4]
    line["salinityAlgo"] = i[8]
    line["densityAlgo"] = i[9]
    line["tempThreshold"] = i[13]
    line["densityThreshold"] = i[17]
    line["tempGradient"] = i[21]
    line["densityGradient"] = i[22]
    out[i[0],i[26]]=line

with open("matOutput.pickle","wb") as f:
    pickle.dump(out,f)

