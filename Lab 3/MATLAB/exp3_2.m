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

fc = 10e3; Ac = 1;
%In-phase carrier
ct1= Ac * cos(2*pi*fc*tt);
%The first message signal
tau_m = 0.001;
mt1 = 2*sinc(tt/tau_m);
%DSB-SC1
st1 = mt1.*ct1;
%Quad. carrier
ct2 = Ac * sin(2*pi*fc*tt);
%The second message signal
mt2 = sinc(tt/tau_m).^2;
%DSB-SC2
st2 = mt2.*ct2;
%QAM signal
st = st1+st2;

%% Plot the QAM Signal in Time Domain
qam_time = figure(1);
time_dom = plot(tt, st, 'LineWidth', 2);
xlabel('Time (s)','FontWeight','bold');
ylabel('\boldmath$s(t)$ \textbf{(V)}','interpreter','latex');
title('QAM Signal: Time Domain');
subtitle('$s(t)$', 'interpreter', 'latex');
axis([-10/fc 10/fc min(st) max(st)]);
qam_time.WindowState = 'maximized';
export_dest = "../Report/Figures/qam_time.png";
exportgraphics(qam_time, export_dest);

%% Plot the QAM Signal in Frequency Domain (Magnitude Spectrum)
Sf = fftshift(fft(fftshift(st)))/(2*fmax);
qam_freq = figure(2);
%The amplitude of the spectrum is different from the Fourier transform
%amplitude due to discretization of discrete Fourier transform
freq_dom = plot(freq,abs(Sf),'LineWidth',2);
xlabel('Frequency (Hz)','FontWeight','bold');
ylabel('\boldmath$|S(f)|$', 'interpreter', 'latex');
title('Spectrum of the QAM Signal S(f)');
subtitle('$|S(f)|$','interpreter','latex');
axis([-30e3 30e3 0 200/(2*fmax)])
qam_freq.WindowState = 'maximized';
export_dest = "../Report/Figures/qam_freq.png";
exportgraphics(qam_freq, export_dest);

%@Rx
%QAM demux and demodulation
% function defined at bottom of script
lpf_signal(st, 0, freq, fmax, mt1, mt2, tt, fc, 1.5e3, 3);
lpf_signal(st, pi/4, freq, fmax, mt1, mt2, tt, fc, 1.5e3, 4);

set(0,'DefaultAxesFontSize','remove');

% function for lpf output
function lpf_signal(st, theta, freq, fmax, mt1, mt2, tt, fc, f_cutoff, fig_num)
    %Local oscillator at the receiver perfectly synchronized
    thet = theta;
    lo1 = cos(2*pi*fc*tt + thet);
    lo2 = sin(2*pi*fc*tt + thet); 
    %upper branch of QAM Rx
    stt1 = st .* lo1;
    Sff1 = fftshift(fft(fftshift(stt1)))/(2*fmax);
    %lower branch of QAM Rx
    stt2 = st .* lo2;
    Sff2 = fftshift(fft(fftshift(stt2)))/(2*fmax);

    %Low pass filtering
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
    Mff1 = Sff1 .* Hf;
    mt_hat1 = 2*fmax*fftshift(ifft(fftshift(Mff1)));
    Mff2 = Sff2 .* Hf;
    mt_hat2 = 2*fmax*fftshift(ifft(fftshift(Mff2)));

    fig = figure(fig_num);
    t = tiledlayout(1,2);
    
    nexttile;
    plt1 = plot(tt,mt_hat1,'r',tt,mt1*0.5,'g.','LineWidth',2);
    xlabel('Time (s)','FontWeight','bold');
    ylabel('\boldmath$\hat{m}_1(t)$ \textbf{(V)}','interpreter','latex')
    title('Output of LPF: Time Domain');
    switch theta
        case 0
            subtitle('$\hat{m}_1(t)$, $\theta = 0$','interpreter','latex');
        case pi/4
            subtitle('$\hat{m}_1(t)$, $\theta = \pi/4$','interpreter','latex');
    end
    axis([-0.01 0.01 min(mt1*0.5) max(mt1*0.5)])
    legend('LPF output', 'Message Signal 1');
    
    nexttile;
    plt2 = plot(tt,mt_hat2,'r',tt,mt2*0.5,'g.','LineWidth',2);
    xlabel('Time (s)','FontWeight','bold');
    ylabel('\boldmath$\hat{m}_2(t)$ \textbf{(V)}','interpreter','latex')
    title('Output of LPF: Time Domain');
    switch theta
        case 0
            subtitle('$\hat{m}_2(t)$, $\theta = 0$','interpreter','latex');
        case pi/4
            subtitle('$\hat{m}_2(t)$, $\theta = \pi/4$','interpreter','latex');
    end
    axis([-0.01 0.01 min(mt2*0.5) max(mt2*0.5)])
    legend('LPF output', 'Message Signal 2');
    
    fig.WindowState = 'maximized';
    export_dest = "../Report/Figures/qam_demod_lpf_theta_" + num2str(thet) + ".png";
    exportgraphics(fig, export_dest);
end