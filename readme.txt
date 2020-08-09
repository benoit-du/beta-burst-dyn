%%% 10-07-20        first revision
%%% Benoit Duchet, University of Oxford
*************************************************
This package can be used to 
- efficiently calculate average burst duration profiles from envelope time series. This is done with fast C code.
- infer envelope models from envelope time series using the passage method. Both the drift function and the noise standard deviation are inferred.
For more details, see B. Duchet, at al., 2020, “Average beta burst duration profiles provide a signature of dynamical changes between the ON and OFF medication states in Parkinson’s disease”, bioRxiv, doi:10.1101/2020.04.27.064246.
*************************************************
*** Run the script avgBurstDurTest.m to test average burst duration profile extraction from the synthetic data provided.

*** Run the script driftInferenceTest.m to test drift inference in the synthetic data provided.

*** Reuse the functions and bits of scripts provided for your own applications. Please send any comment or suggestion to benoit.duchet@ndcn.ox.ac.uk
*************************************************
Note: Non-Windows users should mex the c-files prior to running the scripts. To this end,
	- set 'modules' as the current folder in Matlab
	- run the following commands
			mex fwdSim_learntDyn_mex.c
			mex getBurstDurationAmplitude_mex.c

 
			