%This code sets up the time and frequency vectors for all the numerical
%experiments of Lab3
clear 
format long e
tend = 10;
tbeg = -10;
N=100000;
tstep = (tend-tbeg)/N;
sampling_rate = 1/tstep;

%Time window =
tt = tbeg:tstep:tend-tstep;

%load('lab4_num_expt1')
load('lab4_num_expt2')
%load('lab4_num_expt3')

maxlag = 100;
%Autocorrelation of yt
Ry  = xcorr(yt,yt,maxlag);
%tau vector
tau_vec = -(maxlag*tstep):tstep:maxlag*tstep;
%Abs. PSD corresponding to yt
Sy = abs(fftshift(fft(fftshift(Ry))));
%define the frequency vector corresponding to tau_vec
Ntau = length(tau_vec);
%Nyquist sampling rate
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/Ntau;
%Frequency window
freq = fmin:fstep:fmax-fstep;

hmm = abs(fftshift(fft(fftshift(yt))));
plot(tt,hmm);
