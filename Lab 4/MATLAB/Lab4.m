%This code sets up the time and frequency vectors for all the numerical
%experiments of Lab3
clear 
format long e
set(0,'DefaultAxesFontSize',16);

tend = 10;
tbeg = -10;
N=100000;
tstep = (tend-tbeg)/N;
sampling_rate = 1/tstep;

%Time window =
tt = tbeg:tstep:tend-tstep;

% NUMERICAL EXPERIMENT 1

yt1 = load('../lab4_num_expt1');

lag = [100 200 500];

for i = 1:3
   fig = figure(i);
   tiledlayout(1,2);
   
	% maxlag = 100;
   maxlag = lag(i);
	%Autocorrelation of yt
	Ry  = xcorr(yt1.yt,yt1.yt,maxlag);
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

   nexttile;
   plot(tau_vec, Ry);
   title("Autocorrelation R_y");
   subtitle("maxlag = " + num2str(maxlag));
   xlabel("\tau (s)", 'FontWeight', 'bold');
   ylabel("R_x", 'FontWeight', 'bold');
   ylim([min(Ry) max(Ry)]);
   yline(0);

   nexttile;
   plot(freq, Sy);
   title("Absolute Power Density Function");
   subtitle("maxlag = " + num2str(maxlag));
   xlabel("Frequency (Hz)", 'FontWeight', 'bold');
   ylabel("|S_y|", 'FontWeight', 'bold');
   
   fig.WindowState = 'maximized';
   export_dest = "../Report/Figures/exp1_maxlag_" + num2str(maxlag) + ".png"; 
   exportgraphics(fig, export_dest);
end

% NUMERICAL EXPERIMENT 2

yt2 = load('../lab4_num_expt2');

lag = [100 200 2000];

for i = 4:6
   fig = figure(i);
   tiledlayout(1,2);
   
	% maxlag = 100;
   maxlag = lag(i-3);
	%Autocorrelation of yt
	Ry  = xcorr(yt2.yt,yt2.yt,maxlag);
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

   nexttile;
   plot(tau_vec, Ry);
   title("Autocorrelation vs lag");
   subtitle("maxlag = " + num2str(maxlag));
   yline(0);

   nexttile;
   plot(freq, Sy);
   title("Absolute PSD");
   subtitle("maxlag = " + num2str(maxlag));
   
   fig.WindowState = 'maximized';
end

fig = figure(7);
plot(tt, yt2.yt);
xlim([-100*tstep 100*tstep]);
title("Time Domain");
fig.WindowState = 'maximized';

% NUMERICAL EXPERIMENT 3

exp3_vals = load('../lab4_num_expt3');
exp3_xt = exp3_vals.xt;
exp3_yt = exp3_vals.yt;

fig = figure(8);
t = tiledlayout(1,2);

nexttile;
plot(tt, exp3_xt);
title("xt");

nexttile;
plot(tt, exp3_yt);
title("yt");

fig.WindowState = 'maximized';

maxlag = 20000;
R_xy = xcorr(exp3_yt,exp3_xt,maxlag);
tau_vec = -(maxlag*tstep):tstep:maxlag*tstep;
%Abs. PSD corresponding to yt
S_xy = abs(fftshift(fft(fftshift(R_xy))));
%define the frequency vector corresponding to tau_vec
Ntau = length(tau_vec);
%Nyquist sampling rate
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/Ntau;
%Frequency window
freq = fmin:fstep:fmax-fstep;

fig = figure(9);
plot(tau_vec, R_xy);

[m, i] = max(R_xy);
tau_vec(i)

set(0,'DefaultAxesFontSize','remove');