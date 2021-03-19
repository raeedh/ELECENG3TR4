clear
hold off
format long e
set(0,'DefaultAxesFontSize',16);

N = 4096; %No. of FFT samples
sampling_rate = 100.0e3; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;

tmin = -tmax;
tt = tmin:tstep:tmax-tstep;
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;

% message signal m(t)
fm = 1e3;
Am = 1;
mt = Am * cos(2*pi*fm*tt);
% carrier signal c(t)
fc = 10e3;
Ac = 1;
ct =  Ac * cos(2*pi*fc*tt);
% dsb-sc signal s(t)
st = mt.*ct;

%% Plot DSB-SC Signal in Time Domain
dsbsc_time = figure(1);
time_dom = plot(tt, st, 'LineWidth', 2);
hold on
pos_env = plot(tt,mt,'Color','g','LineStyle','--','LineWidth',2);
neg_env = plot(tt,-mt,'Color','r','LineStyle','--','LineWidth',2);
hold off

% plots phase reversals based on message signal
% function defined at bottom of script
phase_Reversal(tt, tstep, mt, 10/fc, -10/fc);

xlabel('Time (s)','FontWeight','bold');
ylabel('\boldmath$s(t)$ \textbf{(V)}','interpreter','latex');
title('DSB-SC Modulated Wave: Time Domain');
subtitle('$s(t)$', 'interpreter', 'latex');
axis([-10/fc 10/fc min(st) max(st)]);
legend('DSB-SC Signal', 'Positive Envelope', 'Negative Envelope', 'Phase Reversal');
dsbsc_time.WindowState = 'maximized';
export_dest = "../Report/Figures/dsbsc_time.png";
exportgraphics(dsbsc_time, export_dest);

%% Plot DSB-SC in Frequency Domain (Magnitude Spectrum)
Sf = fftshift(fft(fftshift(st)))/(2*fmax);
dsbsc_freq = figure(2);
freq_dom = plot(freq,abs(Sf),'LineWidth',2);
xlabel('Frequency (Hz)','FontWeight','bold');
ylabel('\boldmath$|S(f)|$','interpreter','latex');
title('Magnitude Spectrum of the DSB-SC wave');
subtitle('$|S(f)|$', 'interpreter', 'latex');
axis ([-30e3 30e3 0 max(abs(Sf))]);
dsbsc_freq.WindowState = 'maximized';
export_dest = "../Report/Figures/dsbsc_freq.png";
exportgraphics(dsbsc_freq, export_dest);

%% DSB-SC Demodulation
thet = 0;
lo = cos(2*pi*fc*tt + thet);
st1 = st .* lo;

%% Plot Remodulated Signal in Time Domain
dsbsc_remod_time = figure(3);
time_dom = plot(tt, st1, 'LineWidth', 2);
xlabel('Time (s)','FontWeight','bold');
ylabel('\boldmath$\hat{s}(t)$ \textbf{(V)}','interpreter','latex');
title('DSB-SC Remodulated Wave: Time Domain');
subtitle('$\hat{s}(t)$', 'interpreter', 'latex');
axis([-10/fc 10/fc min(st) max(st)]);
dsbsc_remod_time.WindowState = 'maximized';
export_dest = "../Report/Figures/dsbsc_remod_time.png";
exportgraphics(dsbsc_remod_time, export_dest);

%% Plot Remodulated Signal in Frequency Domain (Magnitude Spectrum)
Sf1 = fftshift(fft(fftshift(st1)))/(2*fmax);
dsbsc_remod_freq = figure(4);
freq_dom = plot(freq,abs(Sf1),'LineWidth',2);
xlabel('Frequency (Hz)','FontWeight','bold');
ylabel('\boldmath$|\hat{S}(f)|$','interpreter','latex');
title('Magnitude Spectrum of the DSB-SC Remodulated Signal');
subtitle('$|\hat{S}(f)|$','interpreter','latex');
axis ([-50e3 50e3 0 max(abs(Sf1))]);
dsbsc_remod_freq.WindowState = 'maximized';
export_dest = "../Report/Figures/dsbsc_remod_freq.png";
exportgraphics(dsbsc_remod_freq, export_dest);

% Low pass filtering
% function defined at end
% 2nd argument is theta (0 for i, pi/2 for ii, pi/4 for iii)
lpf_signal(st, 0, freq, fmax, mt, tt, fc, 5);
lpf_signal(st, pi/2, freq, fmax, mt, tt, fc, 6);
lpf_signal(st, pi/4, freq, fmax, mt, tt, fc, 7);

set(0,'DefaultAxesFontSize','remove');

function Hf = ideal_lpf(f_cutoff, freq)
    %ideal low pass filter
    n=1;
    for f = freq
        if abs(f) < f_cutoff
            Hf(n) = 1;
        else
            Hf(n) = 0;
        end
    n=n+1;
    end
end

function lpf_signal(st, theta, freq, fmax, mt, tt, fc, fig_num)
    % Plots LPF Output for Signals w/ Given Theta (3.1c)
    f_cutoff = 1.5e3;
    Hf = ideal_lpf(f_cutoff, freq);
    thet = theta;
    lo = cos(2*pi*fc*tt + thet);
    st1 = st .* lo;
    Sf1 = fftshift(fft(fftshift(st1)))/(2*fmax);
    Mf1 = Sf1 .* Hf;
    mt1 = 2*fmax*fftshift(ifft(fftshift(Mf1)));
    
    fig = figure(fig_num);
    plt = plot(tt,mt1,'r',tt,mt*0.5,'g--','LineWidth',2);
    xlabel('Time (s)','FontWeight','bold');
    ylabel('\boldmath$\hat{m}(t)$\textbf{ (V)}','interpreter','latex');
    title('Output of LPF: Time Domain');
    switch theta
        case 0
            subtitle("$\hat{m}(t)$, $\theta = 0$", 'interpreter', 'latex');
        case pi/2
            subtitle("$\hat{m}(t)$, $\theta = \pi/2$", 'interpreter', 'latex');
        case pi/4
            subtitle("$\hat{m}(t)$, $\theta = \pi/4$", 'interpreter', 'latex');
    end
    axis([-0.01 0.01 min(mt*0.5) max(mt*0.5)])
    legend('LPF output', 'Message Signal');
    fig.WindowState = 'maximized';
    export_dest = "../Report/Figures/dsbsc_remod_lpf_theta_" + num2str(thet) + ".png";
    exportgraphics(fig, export_dest);
end

function phase_Reversal(tt, tstep, mt, plot_max, plot_min)
    % plot phase reversals based on zeros in message signal
    % zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % find zeros
    zci = @(v) find(diff(sign(v)));
    zst = zci(mt); % find zeros of message signal
    zst = tt(zst);
    for x = zst
        if (x > plot_min) && (x < plot_max)
            if (x > 0) 
                xline(x+tstep,'LineWidth',2,'LineStyle',':');
            else
                xline(x,'LineWidth',2,'LineStyle',':');
            end
        end
    end
end