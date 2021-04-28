import numpy as np
import Sampling_module as Sam
from multiprocessing import Process, Pipe

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

#Rejection sampling with experimental errors.
def Samp(ProcNum, SampleSize, RefSeq, MutSeq, RefProbe, MutProbe, CORp, COMp, Const, log_std, Temp_K, Ratio):

    #Multiprocessing of rejection sampling.
    procs = []
    pipe_list = []
    for p in range(ProcNum):
        recv_end, send_end = Pipe(False)
        proc = Process(target=Sam.Samp, args=(int(SampleSize/ProcNum), RefSeq, MutSeq, RefProbe, MutProbe, CORp, COMp, Const, log_std, Temp_K, Ratio, send_end))
        procs.append(proc)
        pipe_list.append(recv_end)
        proc.start()

    for proc in procs:
        proc.join()

    result_list = np.zeros((9, 0))
    for x in pipe_list:
        result_list = np.concatenate((result_list, np.array(x.recv())), axis=1)

    np.savetxt('CRO_pseudo_list_Sampled_'+str(log_std), result_list[0,:]) #List of sampled initial WT concentrations.
    np.savetxt('CMO_pseudo_list_Sampled_'+str(log_std), result_list[1,:]) #List of sampled initial SNV concentrations.
    np.savetxt('CORp_pseudo_list_Sampled_'+str(log_std), result_list[2,:]) #List of sampled WT concentrations in equilibrium state.
    np.savetxt('COMp_pseudo_list_Sampled_'+str(log_std), result_list[3,:]) #List of sampled SNV concentrations in equilibrium state.

    Ratio_list = np.array(result_list[8,:])

    #Count for True (T) and False (F) samples in +- 0.005 ratio range.
    T, F = 0, 0
    for i in range(np.shape(Ratio_list)[0]):
        if Ratio_list[i] > Ratio - 0.005 and Ratio_list[i] <= Ratio + 0.005:
            T += 1
        else:
            F += 1

    return T, F
