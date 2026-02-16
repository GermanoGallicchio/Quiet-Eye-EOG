% This script demonstrates how to use the QE_EOG function 
% to score QE periods from EOG data

clear; close all; clc;

% Step 1: Add the QE_EOG function to MATLAB path
addpath(pwd);

% Step 2: Load the sample code
load('sample code.mat');

% choose either 3A or 3B
% Step 3A: Choose the algorithm parameters (DISPERSION)
algorithmChoice.name = 'dispersion';
algorithmChoice.winlen = 150;         % Median filter length (samples)
algorithmChoice.threshold = 3;   % numeric (e.g., 3) or 'auto' per-trial


% Step 3B: Choose the algorithm parameters
algorithmChoice.name = 'velocity';
algorithmChoice.winlen = 767;         % Savitzky-Golay frame length (samples)
algorithmChoice.polynDeg = 5;         % Savitzky-Golay polynomial degree
algorithmChoice.threshold = 33;   % numeric (e.g., 33) or 'auto' per-trial

% Step 4: Run the function
% Note: set showPlot = true to visualize the results for the first trial
[QEonset, QEoffset] = QE_EOG(y, x_sec, algorithmChoice, true);

% Step 5: Summarize the results
disp('Quiet Eye Analysis Results:');
disp(sprintf('Number of trials: %d', length(QEonset)));
disp(' ');
disp('QE Onset times (seconds):');
disp(QEonset);
disp(' ');
disp('QE Offset times (seconds):');
disp(QEoffset);
disp(' ');
disp('QE Duration (seconds):');
QE_duration = QEoffset - QEonset;
disp(QE_duration);

% Summary statistics
disp(' ');
disp('Summary Statistics:');
disp(sprintf('Mean QE onset: %.3f sec', nanmean(QEonset)));
disp(sprintf('Mean QE offset: %.3f sec', nanmean(QEoffset)));
disp(sprintf('Mean QE duration: %.3f sec', nanmean(QE_duration)));
disp(sprintf('Std QE duration: %.3f sec', nanstd(QE_duration)));
