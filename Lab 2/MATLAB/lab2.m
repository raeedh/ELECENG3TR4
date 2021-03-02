clear
hold off
format long e

%% Generating time vector, frequency vector, carrier signal, message signal, and modulated signal
N = 2^16; %No. of FFT samples
global sampling_rate % global variable so we can use in functions
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
global fm; % global variable so you can use in functions
fm = 1e3;
Tm = 0.0005;
mt = -2*sinc(tt/Tm);

%% Plotting message signal (Q1)
message_signal = figure(1);
tlayout = tiledlayout(2,1);

% time domain
nexttile;
time_dom = plot(tt, mt, 'LineWidth', 2);
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
freq_dom = plot(freq, abs_Mf, 'LineWidth', 2);
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

% Function for plotting modulated signal in time and frequency domain can
% be found at the bottom of the script
modulated_fig(st, mt, ka, tt, freq, 2, 0.5, "q2");

%% Demodulation (Q2i, Q2ii, Q2iii)

% Envelope detector function can be found at bottom of script
% Function for plotting output of envelope detector and ouput of DC removal
% can be found at bottom of script

% Plotting output of envelope detector and output of DC removal for time constant RC = 1/fc
RC = 1/fc;
yt = envelope_detector(st, RC, tt, N);
output_fig(yt, mt, ka, tt, 3, 0.5, "1/f_c", "q2i", 0);

% Plotting output of envelope detector and output of DC removal for time constant RC = 10*Tm
RC = 10*Tm;
yt = envelope_detector(st, RC, tt, N);
output_fig(yt, mt, ka, tt, 4, 0.5, "10T_m", "q2ii", 0);

% Plotting output for 2iii after exploration of different RC values, the output is passed through a LPF
RC = 1.1*Tm;
yt = envelope_detector(st, RC, tt, N);
output_fig(yt, mt, ka, tt, 5, 0.5, "1.1T_m", "q2iii", 1);

%% Increase modulation to above 100% (Q3)
%For 200% modulation
ka=2/maxmt;

%AM signal
st = (1+ka*mt).*ct;

modulated_fig(st, mt, ka, tt, freq, 6, 2, "q3_mod");

RC = 0.5*(Tm + (1/fc));
yt = envelope_detector(st, RC, tt, N);
output_fig(yt, mt, ka, tt, 7, 2, "0.5(T_m + 1/f_c)", "q3_demod", 1);


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

function modulated_fig(signal, original, ka, time_vector, freq_vector, fig_num, percent_modulation, file_name)
    % Plotting the modulated signal
    fig = figure(fig_num);
    tlayout = tiledlayout(2,1);

    % time domain
    nexttile;
    time_dom = plot(time_vector, signal, 'b', 'LineWidth', 2);
    % plotting envelope
    hold on
    plot(time_vector, ka*original+1,'Color','r','LineStyle', '--','LineWidth',2);
    plot(time_vector, -ka*original-1,'Color','r','LineStyle', '--','LineWidth',2);
    hold off
    legend('Modulated Signal', 'Envelope');
    tim_dom_ax = gca;
    set(tim_dom_ax,'FontSize',16);
    xlabel('Time (s)','FontWeight','bold','Fontsize',16);
    ylabel('Modulated Signal s(t) (V)','FontWeight','bold','Fontsize',16);
    title_name = "Modulated Signal in Time Domain (" + percent_modulation*100 + "% Modulation)";
    title(title_name);
    axis([-2e-3 2e-3 min(signal) max(signal)]);

    % frequency domain
    Sf1 = fft(fftshift(signal));
    Sf = fftshift(Sf1);

    nexttile;
    freq_dom = plot(freq_vector, abs(Sf), 'LineWidth', 2);
    freq_dom_ax = gca;
    set(freq_dom_ax,'FontSize',16);
    xlabel('Frequency (Hz)','FontWeight','bold','Fontsize',16);
    ylabel('|M(f)|','FontWeight','bold','Fontsize',16);
    title_name = "Magnitude Spectrum of the Modulated Signal (" + percent_modulation*100 + "% Modulation)";
    title(title_name);
    axis([-25e3 25e3 0 max(abs(Sf))]);

    fig.WindowState = 'maximized';
    export_dest = "../Report/Figures/" + file_name + ".png";
    exportgraphics(fig, export_dest);
end

function yt = envelope_detector(signal, time_const, time_vector, N)
    tole = 0.1;
    yt = zeros(1,N);
    yt(1) = signal(1);
    n=1;
    for t=time_vector
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

function output_fig(signal, original, ka, time_vector, fig_num, percent_modulation, RC_name, file_name, lpf) % lpf = 0 for no lpf, 1 for lpf
    global sampling_rate;
	global fm;

    fig = figure(fig_num);
    switch lpf
        case 0
            tlayout = tiledlayout(1,2);
        case 1
            tlayout = tiledlayout(1,3);
    end
    title_name = "Output signals for " + percent_modulation*100 + "% Modulation and R_LC = " + RC_name;
    title(tlayout, title_name,'FontWeight','bold','Fontsize',24);

    % output of envelope detector
    nexttile;
    envelope_det = plot(time_vector, signal,'LineWidth',2);
    envelope_det_ax = gca;
    set(envelope_det_ax,'FontSize',16);
    xlabel('Time (s)','FontWeight','bold','Fontsize',16);
    ylabel('y(t) (V)','FontWeight','bold','Fontsize',16);
    title('After the envelope detector');
    axis([-2e-3 2e-3 0 max(signal)]);

    % dc removal and division by ka
    yt1 = (signal - 1) / ka;

    nexttile;
    output_signal = plot(time_vector,yt1,'r',time_vector,original,'k','LineWidth',2);
    legend('after DC removal','message signal');
    output_signal_ax = gca;
    set(output_signal_ax,'FontSize',16);
    xlabel('Time (s)','FontWeight','bold','Fontsize',16);
    ylabel('y1(t) (V)','FontWeight','bold','Fontsize',16);
    title('After the DC removal');
    axis([-2e-3 2e-3 min(original) max(original)]);
    
    if (lpf == 1)
        nexttile;
        filtered_signal = plot(time_vector,lowpass(yt1, 1.1*fm, sampling_rate),'r',time_vector,original,'k','LineWidth',2);
        legend('after LPF','message signal');
        filtered_signal_ax = gca;
        set(filtered_signal_ax,'FontSize',16);
        xlabel('Time (s)','FontWeight','bold','Fontsize',16);
        ylabel('m1(t) (V)','FontWeight','bold','Fontsize',16);
        title('After the LPF');
        axis([-2e-3 2e-3 min(original) max(original)]);
    end

    fig.WindowState = 'maximized';
    export_dest = "../Report/Figures/" + file_name + ".png";
    exportgraphics(fig, export_dest);
end
