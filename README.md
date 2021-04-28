# MismatchIntroducedProbe
Files and codes list and discription.

1. Run_Mq(x).py
Code for generation of linear regression parameters of q(x).
q(x) parameters should be generated before the actual rejection sampling.

2. Run_Sampling.py
Reject sampling execution code.

3. Sampling_Multi.py
Code for preparation and collection of accepted samples vai python multiprocessing procedure.

4. Sampling_module.py
Code for rejection sampling. Imported in Sampling_Multi.py code and executed for the sampling.

5. nupack_wrapper.py
Code for execution of nupack software for calculation of equilibrium state.

6. MutantSequeces (Directory)
Sequence of WT and SNVs included in.
