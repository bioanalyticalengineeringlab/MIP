import numpy as np
import Hybridization_NUPACK as HyN
from Bio.Seq import Seq as Sq
import glob
import os

fileList = glob.glob('/tmp/tmp*')
for filePath in fileList:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)

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


#Rejection sampling
import Sampling_Multi as Sam

log_std = -10 #log scaled standard deviation from experimental error.

SampleSize = 1000
ProcNum = 50

Ratio = CMO / (CMO + CRO)

T, F = Sam.Samp(ProcNum, SampleSize, RefSeq, MutSeq, RefProbe, MutProbe, CORp, COMp, Const, log_std, Temp_K, Ratio)
