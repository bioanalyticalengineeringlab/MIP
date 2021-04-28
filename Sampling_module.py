import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import Hybridization_NUPACK as HyN

#Import q(x) parameters.
def Mq_Rx(x, c):
    para_Rx = np.loadtxt('popt_CORp_pseudo')
    a = para_Rx[0]
    b = para_Rx[1]

    if c == 0:
        return ((a * x) + b)
    else:
        return ((a * x) + b) * c

def Mq_Mx(x, c):
    para_Mx = np.loadtxt('popt_COMp_pseudo')
    a = para_Mx[0]
    b = para_Mx[1]

    if c == 0:
        return ((a * x) + b)
    else:
        return ((a * x) + b) * c

def Samp(SampleSize, RefSeq, MutSeq, RefProbe, MutProbe, CORp, COMp, Const, log_std, Temp_K, Ratio, send_end):
    np.random.seed(None)
    std = math.pow(10, log_std)

    #Adjustment of q(x) parameters for making q(x) always higher than p(x).
    Rx_adjust = (stats.norm(CORp, std).pdf(CORp) * 3 / 2) / Mq_Rx(CORp, 0)
    Mx_adjust = (stats.norm(COMp, std).pdf(COMp) * 3 / 2) / Mq_Mx(COMp, 0)

    CRO_pseudo_list_Sampled = []
    CMO_pseudo_list_Sampled = []
    CORp_pseudo_list = []
    COMp_pseudo_list = []
    px_R_list = []
    px_M_list = []
    Mq_R_list = []
    Mq_M_list = []

    done = False
    while (not done):
        #Random sample generation.
        CRO_pseudo_list = np.random.uniform(0, Const, 1)
        CMO_pseudo_list = Const - CRO_pseudo_list

        CRO_pseudo, CMO_pseudo = CRO_pseudo_list[0], CMO_pseudo_list[0]

        # Calculation of equilibrium states of hybridization.
        RefC_RefProbe, MutC_RefProbe, RefG_RefProbe, MutG_RefProbe = HyN.Hyb(RefSeq, MutSeq, RefProbe, CRO_pseudo,
                                                                             CMO_pseudo, Temp_K)
        RefC_MutProbe, MutC_MutProbe, RefG_MutProbe, MutG_MutProbe = HyN.Hyb(RefSeq, MutSeq, MutProbe, CRO_pseudo,
                                                                             CMO_pseudo, Temp_K)

        CORp_pseudo = RefC_RefProbe + MutC_RefProbe
        COMp_pseudo = RefC_MutProbe + MutC_MutProbe

        px_R, px_M = stats.norm(CORp, std).pdf(CORp_pseudo), stats.norm(COMp, std).pdf(COMp_pseudo)
        Mq_R, Mq_M = Mq_Rx(CORp_pseudo, Rx_adjust), Mq_Mx(COMp_pseudo, Mx_adjust)

        Threshold_R = np.random.uniform(0, Mq_R, 1)[0]
        Threshold_M = np.random.uniform(0, Mq_M, 1)[0]

        #Acceptance determination.
        if px_R >= Threshold_R and px_M >= Threshold_M:
            CRO_pseudo_list_Sampled.append(CRO_pseudo)
            CMO_pseudo_list_Sampled.append(CMO_pseudo)

            CORp_pseudo_list.append(CORp_pseudo)
            COMp_pseudo_list.append(COMp_pseudo)

            px_R_list.append(px_R)
            px_M_list.append(px_M)
            Mq_R_list.append(Mq_R)
            Mq_M_list.append(Mq_M)

            print('Sampled', len(px_R_list))
            if len(px_R_list) >= SampleSize:
                done = True

    plt.subplot(4, 1, 1)
    plt.plot(CORp_pseudo_list, px_R_list, 'ro')

    plt.subplot(4, 1, 2)
    plt.plot(CORp_pseudo_list, Mq_R_list, 'bo')

    plt.subplot(4, 1, 3)
    plt.hist(px_R_list, 100)

    plt.subplot(4, 1, 4)
    plt.hist(Mq_R_list, 100)


    Ratio_list = np.array(CMO_pseudo_list_Sampled) / (np.array(CMO_pseudo_list_Sampled) + np.array(CRO_pseudo_list_Sampled))

    send_end.send(np.array(
        [CRO_pseudo_list_Sampled, CMO_pseudo_list_Sampled, CORp_pseudo_list, COMp_pseudo_list, px_R_list,
         px_M_list, Mq_R_list, Mq_M_list, Ratio_list]))
