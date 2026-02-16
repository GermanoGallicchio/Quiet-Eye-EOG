# QE_EOG
MATLAB code for scoring Quiet Eye durations from EOG data

how to apply this function to the sample data
1. add the QE_EOG function to MATLAB path, for example using addpath()
2. load the sample code
3. choose the algorithm parameters, for example: 
   algorithmChoice.name = 'velocity';
   algorithmChoice.winlen = 767;        % Savitzky-Golay frame length (samples)
   algorithmChoice.polynDeg = 5;        % Savitzky-Golay polynomial degree
   algorithmChoice.threshold = 33;  % numeric (e.g., 33) or 'auto' per-trial
4. run the function: QE_EOG(y, x_sec, algorithmChoice, true)

Notes
- Inputs can be single-trial (vector) or multi-trial (samples × trials).
- `threshold` accepts a scalar, a per-trial vector, or `'auto'` to estimate a per-trial threshold (stable-window robust variance → 3× std).
- `showPlot = true` plots only the first trial when multi-trial input is provided.
- `x_sec` must span negative to positive time and include exactly one 0 (movement initiation reference).

## Cite as
Gallicchio, G. (2023). Quiet-Eye-EOG. _Zenodo_. https://doi.org/10.5281/zenodo.8411092
