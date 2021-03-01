clear
hold off
format long e

%% Generating time vector, frequency vector, carrier signal, message signal, and modulated signal
N = 2^16; %No. of FFT samples
sampling_rate = 40e4; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;
tmin = -tmax;
tt = tmin:tstep:tmax-tstep;
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;

%carrier
fc=20e3;
Ac = 1;
ct=Ac*cos(2*pi*fc*tt);

%message signal
fm = 1e3;
Tm = 0.0005;
mt = -2*sinc(tt/Tm);

%% Plotting message signal (Q1)
message_signal = figure(1);
tlayout = tiledlayout(2,1);

% time domain
nexttile;
time_dom = plot(tt, mt);
set(time_dom,'LineWidth',2);
tim_dom_ax = gca;
set(tim_dom_ax,'FontSize',16);
xlabel('Time (s)','FontWeight','bold','Fontsize',16);
ylabel('Message m(t) (V)','FontWeight','bold','Fontsize',16);
title('Message Signal in Time Domain');
axis([-2e-3 2e-3 min(mt) max(mt)]);

% frequency domain
Mf1 = fft(fftshift(mt));
Mf = fftshift(Mf1);
abs_Mf = abs(Mf);

nexttile;
freq_dom = plot(freq, abs_Mf);
set(freq_dom,'LineWidth',2);
freq_dom_ax = gca;
set(freq_dom_ax,'FontSize',16);
xlabel('Frequency (Hz)','FontWeight','bold','Fontsize',16);
ylabel('|M(f)|','FontWeight','bold','Fontsize',16);
title('Magnitude Spectrum of the Message Signal');
axis ([-5e3 5e3 0 max(abs(Mf))]);

message_signal.WindowState = 'maximized';
exportgraphics(message_signal,'../Report/Figures/q1.png');

% Measure highest frequency component from positive side of spectrum
[max_val, max_freq_index] = max(abs_Mf((length(abs_Mf)/2+1):end));
highest_freq = freq(max_freq_index + length(freq)/2);

%% Modulated Signal (Q2)

%max of absolute of m(t)
maxmt = max(abs(mt));
%For 50% modulation
ka=0.5/maxmt;

%AM signal
st = (1+ka*mt).*ct;

% Plotting the modulated signal
modulated_signal = figure(2);
tlayout = tiledlayout(2,1);

% time domain
nexttile;
time_dom = plot(tt, st);
set(time_dom,'LineWidth',2);
tim_dom_ax = gca;
set(tim_dom_ax,'FontSize',16);
xlabel('Time (s)','FontWeight','bold','Fontsize',16);
ylabel('Modulated Signal s(t) (V)','FontWeight','bold','Fontsize',16);
title('Modulated Signal in Time Domain');
axis([-2e-3 2e-3 min(st) max(st)]);

% frequency domain
Sf1 = fft(fftshift(st));
Sf = fftshift(Sf1);

nexttile;
freq_dom = plot(freq, abs(Sf));
set(freq_dom,'LineWidth',2);
freq_dom_ax = gca;
set(freq_dom_ax,'FontSize',16);
xlabel('Frequency (Hz)','FontWeight','bold','Fontsize',16);
ylabel('|M(f)|','FontWeight','bold','Fontsize',16);
title('Magnitude Spectrum of the Modulated Signal');
axis([-25e3 25e3 0 max(abs(Sf))]);

modulated_signal.WindowState = 'maximized';
exportgraphics(message_signal,'../Report/Figures/q2.png');

%% Demodulation

% Envelope detector function can be found at bottom of script

% Plotting output of envelope detector and output of DC removal for time constant RC = 1/fc
RC = 1/fc;
yt = envelope_detector(st, RC, tt, N);
plot_2i = figure(3);
tlayout = tiledlayout(1,2);
title(tlayout,'Output signals for R_LC = 1/f_c','FontWeight','bold','Fontsize',24);

% output of envelope detector
nexttile;
envelope_det = plot(tt, yt);
set(envelope_det,'LineWidth',2);
envelope_det_ax = gca;
set(envelope_det_ax,'FontSize',16);
xlabel('Time (s)','FontWeight','bold','Fontsize',16);
ylabel('y(t) (V)','FontWeight','bold','Fontsize',16);
title('After the envelope detector');
axis([-2e-3 2e-3 0 max(yt)]);

% dc removal and division by ka
yt1 = (yt - 1) / ka;

nexttile;
output_signal = plot(tt,yt1,'r',tt,mt,'k');
set(output_signal,'LineWidth',2);
legend('after DC removal','message signal')
output_signal_ax = gca;
set(output_signal_ax,'FontSize',16);
xlabel('Time (s)','FontWeight','bold','Fontsize',16);
ylabel('y1(t) (V)','FontWeight','bold','Fontsize',16);
title('After the DC removal');
axis([-2e-3 2e-3 min(mt) max(mt)]);

plot_2i.WindowState = 'maximized';
exportgraphics(plot_2i,'../Report/Figures/q2i.png');

% Plotting output of envelope detector and output of DC removal for time constant RC = 10*Tm
RC = 10*Tm;
yt = envelope_detector(st, RC, tt, N);
plot_2ii = figure(4);
tlayout = tiledlayout(1,2);
title(tlayout,'Output signals for R_LC = 10T_m','FontWeight','bold','Fontsize',24);

% output of envelope detector
nexttile;
envelope_det = plot(tt, yt);
set(envelope_det,'LineWidth',2);
envelope_det_ax = gca;
set(envelope_det_ax,'FontSize',16);
xlabel('Time (s)','FontWeight','bold','Fontsize',16);
ylabel('y(t) (V)','FontWeight','bold','Fontsize',16);
title('After the envelope detector');
axis([-2e-3 2e-3 0 max(yt)]);

% dc removal and division by ka
yt1 = (yt - 1) / ka;

nexttile;
output_signal = plot(tt,yt1,'r',tt,mt,'k');
set(output_signal,'LineWidth',2);
legend('after DC removal','message signal')
output_signal_ax = gca;
set(output_signal_ax,'FontSize',16);
xlabel('Time (s)','FontWeight','bold','Fontsize',16);
ylabel('y1(t) (V)','FontWeight','bold','Fontsize',16);
title('After the DC removal');
axis([-2e-3 2e-3 min(mt) max(mt)]);

plot_2ii.WindowState = 'maximized';
exportgraphics(plot_2ii,'../Report/Figures/q2ii.png');

% figure(1)
% Hp1 = plot(tt,ct);
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% Hx=xlabel('Time (sec) ');
% set(Hx,'FontWeight','bold','Fontsize',16)
% Hx=ylabel('Carrier c(t)  (Volt)');
% set(Hx,'FontWeight','bold','Fontsize',16)
% title('Carrier : Time domain');
% axis([-1e-3 1e-3 -1.1 1.1])

% %pause
% %Low pass filter to remove the ripple
% %choose the cutoff frequency of the filter to be slightly higher than
% %the highest freq of the message signal
% % f0 = 1.1*fm;
% % mt1 = rect_filt(yt1,freq,f0);
% % figure(8)
% % Hp1 = plot(tt,mt1);
% % set(Hp1,'LineWidth',2)
% % Ha = gca;
% % set(Ha,'Fontsize',16)
% % Hx=xlabel('Time (sec) ');
% % set(Hx,'FontWeight','bold','Fontsize',16)
% % Hx=ylabel('m1(t)  (Volt)');
% % set(Hx,'FontWeight','bold','Fontsize',16)
% % title('After the low pass filter');
% % axis([-2e-3 2e-3 min(mt1) max(mt1)])

function yt = envelope_detector(signal, time_const, tt, N)
    tole = 0.1;
    yt = zeros(1,N);
    yt(1) = signal(1);
    n=1;
    for t=tt
        if(n > 1)
         if(signal(n) > yt(n-1))
             yt(n) = signal(n);
         else
             if((yt(n-1)-signal(n)) < tole)
               yt0 = yt(n-1);
               yt(n) = yt0;
               %time when C starts discharging
               tc = t;
             else
               yt(n) = yt0*exp(-(t-tc)/time_const);
             end
         end
        end
        n=n+1;
    end
    yt(1)=yt(2);
end