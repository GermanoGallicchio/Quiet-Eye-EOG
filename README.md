# QE_EOG
MATLAB code for scoring Quiet Eye durations from EOG data

how to apply this function to the sample data
1. add the QE_EOG function to MATLAB path, for example using addpath()
2. load the sample code
3. choose the algorithm parameters, for example: 
   algorithmChoice.name = 'velocity';
   algorithmChoice.winlen = 767;
   algorithmChoice.polynDeg = 5;
   algorithmChoice.threshold = 33;
4. run the function: QE_EOG(y,x_sec,algorithmChoice,"yes")

## Cite as
https://doi.org/10.5281/zenodo.8411093
