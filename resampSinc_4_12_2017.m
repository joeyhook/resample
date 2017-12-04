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
clear, close

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% User defined parameters:
%   step:   Less than zero time stretches and shifts pitch down,
%           greater than zero time compresses and shifts pitch up.
%           Integers are half steps.
%   f1:     Frequency before pitchshift
%   Tf:     Length of block in seconds
%   Fs_out: Original sample rate that pitchshifted signal is scaled back to
%   sw:     Sinc width - number of zero crossings in sinc function
%           on each side of target sample
%   alpha:  Kaiser window parameter
%   k_sinc:     resolution of table
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step = -7;      % any real number, integers are half steps
f1 = 1.3;       % > 0
Tf = 0.76;      % > 0
Fs_out = 81;    % 44100 is standard, ~50 is good for viewing plot
sw = 7;         % integer, >= 5 is best
alpha = 22;     % 22 seems to work OK
k_sinc = 0.001; % 1/(10^x)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Dependent parameters
%   Q:     Pitchshift ratio (also Fs_in/Fs_out) from # steps
%   k_out: Time-step size of resampled output on seconds
%   N_out: Number of samples in  block
%   t1:    Time vector for plotting before pitchshift and after
%          resample
%   Fs_in: Sample rate after pitchshift / before resampling
%   k_in:  time-step size after pitchshift, before resampling
%   N_in:  Number of samples in block
%   t2:    Time vector for plotting preshift and resampout
%   f2:    Shifted frequency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
semitone = nthroot(2,12);
Q = semitone^step;

k_out  = 1/Fs_out;
N_out  = floor(Tf/k_out);
t1   = 0:k_out:Tf-k_out;

Fs_in = Q*Fs_out;
k_in  = 1/Fs_in;
N_in  = floor(Tf/k_in);
t2   = 0:k_in:Tf-k_in;

f2 = Q*f1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Create original signal and virtual shifted signal vectors
% Column 1 is audio signal data and column 2 is time data
% Output vector initialized with zeros for audio data and
%   time data of preshifted signal   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
preshift = zeros(N_out,2);
for n=1:N_out
    preshift(n,1) = sin(2*pi*f1*(n-1)*k_out);
    preshift(n,2) = (n-1)*k_out;
end

resampin = zeros(N_in,2);
for n=1:N_in
    resampin(n,1) = sin(2*pi*f2*(n-1)*k_in);
    resampin(n,2) = (n-1)*k_in;
end

resampout = zeros(N_out,2);
resampout(:,2) = preshift(:,2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute table of Kaiser windowed sinc finction
% This will only need to be done once per block, and only if 
% Fs_out, sw, or alpha changes
%   k_sinc:     resolution of table (defined in user variables)
%   tabmult:    multiplier to get an integer to call table
%   sw:         Sinc width - number of zero crossings in sinc function
%               on each side of target sample (defined in user variables)
%   t_sinc:     time vector of sinc table
%   t_kais:     time vector of kaiser window
%   N_sinc:     size of sinc table
%   N_kais:     size of Kaiser table
%   piFst:      pre-calculation for input to sinc
%   sinctab:    compute table (not windowed)
%   alpha:      Kaiser window parameter (defined in user variables)
%   kais:       Kaiser window table
%   sinctabwin: wondowed sinc table
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tablestart=tic; % time table creation
tabmult = 1/k_sinc;
t_sinc = 0:k_sinc:sw-k_sinc; % sinc is even, so we only need half
t_kais = -sw:k_sinc:sw-k_sinc; % but MATLAB won't do half of Kaiser
N_sinc = length(t_sinc);
N_kais = length(t_kais);
piFst = (pi*Fs_out).*t_sinc;
oneOverpiFst = 1./piFst;
sinctab = (sin(pi*t_sinc)./(pi*t_sinc))';
sinctab(1) = 1; % was NaN
kais = kaiser(N_kais,alpha);
sinctabwin(1:N_sinc) = sinctab(1:N_sinc).*kais(0.5*N_kais+1:N_kais);
sinctabwin(N_sinc+1)=0; % one larger than sinctab to not go out of bounds
toc(tablestart)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  PLOT first 2 waveforms
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(3,2,1), stem(t1,preshift(:,1)), axis([0,Tf,-1,1])
title('before pitchshift')
subplot(3,2,3), stem(t2,resampin(:,1)), axis([0,Tf,-1,1]),
title('after pitchshift/resample input')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Main loop
%   sincin:     stores values going into sinc table
%   sincoutwin: stores values coming out of table   
%   n_0 & t_0: sample in output to calculate and its time
%   n_inm1 & n_inp1: samples just before and after n_0 from input
%   t_inm1 & t_inp1: time in seconds at n_inm1 & n_inp1
%   eta_L & eta_R:  time between n_0 and next sample of input
%       used to find spot on sinc table
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sincin = zeros(2*sw,1);
sincoutwin = zeros(2*sw,1);
resampout(1) = resampin(1); % first value in block will be the same
mainloopstart = tic; % time main loop
for n_0=2:N_out-1 % -1 is cheap fix for issue 3.  
    t_0 = resampout(n_0,2);
    n_inm1 = floor(t_0*Fs_in)+1;
    n_inp1 = n_inm1 + 1;
    t_inm1 = resampin(n_inm1,2);
    t_inp1 = resampin(n_inp1,2);
    eta_L = t_0 - t_inm1;
    eta_R = t_inp1 - t_0;
    % where to pluck values out of table for left of n_0
    for m=1:sw
        sincin(sw+1-m) = round(tabmult*(Fs_in*eta_L+m-1))+1;
    end
    % where to pluck values out of table for right of n_0
    for m=(sw+1):2*sw
        sincin(m) = round(tabmult*(Fs_in*eta_R+m-sw-1))+1;
    end
    
    % get values from table
    sincoutwin = sinctabwin(sincin);

    % calculate x(n_0)
    for m=1:sw
        indexL = n_inm1-sw+m;
        if indexL < 1
            indexL = 1;
        end
        resampout(n_0) = resampout(n_0) + resampin(indexL)*sincoutwin(m);
    end
    for m=(sw+1):2*sw
        indexR = n_inm1-sw+m;
        if indexR > N_in
            indexR = N_in;
        end
        resampout(n_0) = resampout(n_0) + resampin(indexR)*sincoutwin(m);    
    end
end
toc(mainloopstart)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  PLOT the rest
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(3,2,5), stem(t1,resampout(:,1)), title('resample output')
axis([0,Tf,-1,1]), 

subplot(3,2,4), plot(t_kais,kais), xlim([0,sw]), title('kaiser window')
t_sinc = 0:k_sinc:sw;
subplot(3,2,6), plot(t_sinc,sinctabwin), xlim([0,sw]), title('windowed sinc')
hold on
stem(sincin*k_sinc,sincoutwin)
hold off

