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
	fmax = sampling_rate/2; 
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
    ylabel("R_y", 'FontWeight', 'bold');
    ylim([min(Ry) max(Ry)]);
    yline(0);

    nexttile;
    plot(freq, Sy);
    title("Absolute Power Spectral Density");
    subtitle("maxlag = " + num2str(maxlag));
    xlabel("Frequency (Hz)", 'FontWeight', 'bold');
    ylabel("|S_y|", 'FontWeight', 'bold');
    
    fig.WindowState = 'maximized';
    export_dest = "../Report/Figures/exp1_maxlag_" + num2str(maxlag) + ".png";
    exportgraphics(fig, export_dest);
end

% NUMERICAL EXPERIMENT 2

yt2 = load('../lab4_num_expt2');

lag = [100 200 20000];

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
	fmax = sampling_rate/2; 
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
    ylabel("R_y", 'FontWeight', 'bold');
    ylim([min(Ry) max(Ry)]);
    yline(0);

    nexttile;
    plot(freq, Sy);
    title("Absolute Power Spectral Density");
    subtitle("maxlag = " + num2str(maxlag));
    xlabel("Frequency (Hz)", 'FontWeight', 'bold');
    ylabel("|S_y|", 'FontWeight', 'bold');

    f = Sy(length(Sy)/2:end);
    [m, i] = max(f);
    freq(i+length(f)-1)

    fig.WindowState = 'maximized';
    export_dest = "../Report/Figures/exp2_maxlag_" + num2str(maxlag) + ".png";
    exportgraphics(fig, export_dest);
end

fig = figure(7);
plot(tt, yt2.yt);
xlim([-100*tstep 100*tstep]);
title("Signal yt(n): Time Domain");
subtitle("$y(t)$", 'interpreter', 'latex')
xlabel("Time (s)", 'FontWeight', 'bold');
ylabel("y(t)", 'FontWeight', 'bold');

fig.WindowState = 'maximized';
export_dest = "../Report/Figures/exp2_time.png";
exportgraphics(fig, export_dest);

fig = figure(8);

Nyt = length(yt2.yt);
%Nyquist sampling rate
fmax = sampling_rate/2; 
fmax = sampling_rate/2; 
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/Nyt;
%Frequency window
freq = fmin:fstep:fmax-fstep;

plot(freq, abs(fftshift(fft(fftshift(yt2.yt)))));
title("Signal yt(n): Magnitude Spectrum");
xlabel("Frequency (Hz)", 'FontWeight', 'bold');
ylabel("Magnitude", 'FontWeight', 'bold');

fig.WindowState = 'maximized';
export_dest = "../Report/Figures/exp2_fft.png";
exportgraphics(fig, export_dest);

% NUMERICAL EXPERIMENT 3

exp3_vals = load('../lab4_num_expt3');
exp3_xt = exp3_vals.xt;
exp3_yt = exp3_vals.yt;

fig = figure(9);
t = tiledlayout(1,2);

nexttile;
plot(tt, exp3_xt);
title("Signal x: Time Domain");
subtitle("$x(t)$", 'interpreter', 'latex');
xlabel("Time (s)", 'FontWeight', 'bold');
ylabel("x(t)", 'FontWeight', 'bold');

nexttile;
plot(tt, exp3_yt);
title("Signal y: Time Domain");
subtitle("$y(t)$", 'interpreter', 'latex');
xlabel("Time (s)", 'FontWeight', 'bold');
ylabel("y(t)", 'FontWeight', 'bold');

fig.WindowState = 'maximized';
export_dest = "../Report/Figures/exp3_time.png";
exportgraphics(fig, export_dest);

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

fig = figure(10);
plot(tau_vec, R_xy);
title("Autocorrelation R_{xy}");
subtitle("maxlag = " + num2str(maxlag));
xlabel("\tau (s)", 'FontWeight', 'bold');
ylabel("R_{xy}", 'FontWeight', 'bold');

fig.WindowState = 'maximized';
export_dest = "../Report/Figures/exp3_crosscorr.png";
exportgraphics(fig, export_dest);

[m, i] = max(R_xy);
tau_vec(i)

set(0,'DefaultAxesFontSize','remove');