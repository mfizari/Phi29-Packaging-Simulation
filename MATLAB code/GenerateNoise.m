function [noise_output] = GenerateNoise(cfold, windowsize_input, fname, units)

%This function generates noise in a window of windowsize_input-length (s)
%"Noise" is defined as deviations from the mean of a windowsize_noise window 

%Output as given as bp/kbp with a sr of 1e3

%   Inputs :
%       windowsize : the window size of input data to which to add noise (s) 
%       windowsize_noise : the window size used for generating the noise dist (s)
%       units : 'bp' for basepairs, anything else for kbp

%Determine units to use for noise output
if strcmp(units, 'bp')
    k_units = 1;
else
    k_units = 1e-3;
end

%Navigate to code folder
cd(cfold) 

%Load distribution
load(fname, 'pd'); %this loads "pd" as the kernel distribution for noise

%Generate noise and output!
noise_output = k_units*random(pd, 1e3*windowsize_input,1);
