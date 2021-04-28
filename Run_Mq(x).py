import numpy as np
import Hybridization_NUPACK as HyN
from Bio.Seq import Seq as Sq

#Theoretical hybridization states calculation.
Temp_K = 298

CodeSet = np.loadtxt('./MutantSequences/L858RsetCode')
SeqSet = open('./MutantSequences/L858Rset', 'r')
Seqs = SeqSet.readlines()

Const = 1e-9 #Total concentration of WT and SNV
CRO = Const
CMO = CRO * 0.01 #Concentration of initial SNV
CRO = Const - CMO #Concentration of initial WT

RefCode, MutCode = CodeSet[0, :], CodeSet[1, :]
RefSeq, MutSeq = Seqs[0], Seqs[1] #Sequence of WT and SNV

RefProbe = str(Sq.reverse_complement(Sq(RefSeq))).strip('\n')
MutProbe = str(Sq.reverse_complement(Sq(MutSeq))).strip('\n')

#Calculation of equilibrium states of hybridization.
RefC_RefProbe, MutC_RefProbe, RefG_RefProbe, MutG_RefProbe = HyN.Hyb(RefSeq, MutSeq, RefProbe, CRO, CMO, Temp_K)
RefC_MutProbe, MutC_MutProbe, RefG_MutProbe, MutG_MutProbe = HyN.Hyb(RefSeq, MutSeq, MutProbe, CRO, CMO, Temp_K)

CORp = RefC_RefProbe + MutC_RefProbe
COMp = RefC_MutProbe + MutC_MutProbe

print(CORp, COMp)

#Generation of samples for q(x) regression.
import matplotlib.pyplot as plt
import os
import glob

SampleSize = 10000

CRO_pseudo_list = np.random.uniform(0,Const,SampleSize)
CMO_pseudo_list = Const - CRO_pseudo_list


CORp_pseudo_list = []
COMp_pseudo_list = []
for i in range(np.shape(CRO_pseudo_list)[0]):
    print(i)
    CRO_pseudo, CMO_pseudo = CRO_pseudo_list[i], CMO_pseudo_list[i]

    RefC_RefProbe, MutC_RefProbe, RefG_RefProbe, MutG_RefProbe = HyN.Hyb(RefSeq, MutSeq, RefProbe, CRO_pseudo, CMO_pseudo, Temp_K)
    RefC_MutProbe, MutC_MutProbe, RefG_MutProbe, MutG_MutProbe = HyN.Hyb(RefSeq, MutSeq, MutProbe, CRO_pseudo, CMO_pseudo, Temp_K)

    CORp_pseudo = RefC_RefProbe + MutC_RefProbe
    COMp_pseudo = RefC_MutProbe + MutC_MutProbe

    CORp_pseudo_list.append(CORp_pseudo)
    COMp_pseudo_list.append(COMp_pseudo)

print(CORp_pseudo_list)
print(COMp_pseudo_list)
count, bins, ignored = plt.hist(CORp_pseudo_list, bins=int(SampleSize *0.3))
count2, bins2, ignored2 = plt.hist(COMp_pseudo_list, bins=int(SampleSize *0.3))

np.savetxt('CRO_pseudo', CRO_pseudo_list) #List of sampled initial WT concentration.
np.savetxt('CORp_pseudo', np.array(CORp_pseudo_list)) #List of sampled WT concentration in equilibrium state.
np.savetxt('COMp_pseudo', np.array(COMp_pseudo_list)) #List of sampled SNV concentration in equilibrium state.

#Linear regression of q(x).
from scipy import optimize

def func(x, a, b):
    return (a * x) + b

bin_centers = bins[:-1] + np.diff(bins) / 2
popt_CORp_pseudo, pcov = optimize.curve_fit(func, bin_centers, count)
print(popt_CORp_pseudo)

bin_centers2 = bins2[:-1] + np.diff(bins2) / 2
popt_COMp_pseudo, pcov2 = optimize.curve_fit(func, bin_centers2, count2)
print(popt_COMp_pseudo)

plt.plot(bin_centers, func(bin_centers, *popt_CORp_pseudo), color='red', linewidth=2)
plt.plot(bin_centers2, func(bin_centers2, *popt_COMp_pseudo), color='blue', linewidth=2)


np.savetxt('popt_CORp_pseudo', popt_CORp_pseudo) #Parameters in linear regression of WT.
np.savetxt('popt_COMp_pseudo', popt_COMp_pseudo) #Parameters in linear regression of SNV.

#Remove Temporary files.

fileList = glob.glob('/tmp/tmp*')
for filePath in fileList:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)