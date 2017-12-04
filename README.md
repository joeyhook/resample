# resample
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Resampler using bandlimited interpolation (sinc interpolation)
% Joey Hook - 4/12/2017
%
%   Generates a sinusoid of specified frequency, duration, and
%   sampling rate, then generates another at Q times freq and SR,
%   then resamples the 2nd sinusoid back to the 1st SR.
%   The 3 sinsoids are plotted, along with the Kaiser window and 
%   sinc table used in the interpolation.
%
%   Issues:
%      1. Need ideal value for Kaiser parameter alpha to match
%       output amplitude to input amplitude
%      2. Calculate Kaiser window from scratch for use outside
%       MATLAB and more efficiency
%      3. Problems at end of block when input sample times exceed
%       existing samples. Will be addressed in real-time
%       implementation.   
%      4. indexL and indexR do not go outside block as they should
%       This will also be left to fix in implementation.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
