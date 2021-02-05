clc
clear all

%% Input plots
figure(1);
input_plots = tiledlayout(2,2);
title(input_plots, 'Input Plots', 'FontSize', 20);

%% Square wave generator (input)
f0 = 10000; % fundamental freq of input square wave
T0 = 1/f0;  % period
tstep = 0.001*T0;
no_sample = 6*T0/tstep + 1;
tt = -3*T0:tstep:3*T0;
fs = 1/tstep;   % sampling frequency

input = square(tt*2*pi*f0,50); % input square wave

% plot square wave
nexttile([1 2]);
Hp1 = plot(tt,input);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
xlim([tt(1) tt(length(tt))]);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Input - Time Domain')

%% Fourier series representation of signal (Amplitude Spectrum)
N = 100; % number of harmonics
nvec = -N:N;
c_in = zeros(size(nvec));
for n = nvec
   m = n + N + 1;
   if (mod(n,2))
       % c_in(m) = sinc(n/2);
       c_in(m) = (2/(n*pi))*sin((n*pi)/2);
   else
       c_in(m) = 0.0;
   end
end
f = nvec*f0;
nexttile;
Hp1=stem(f,abs(c_in));
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Magnitude Spectrum of Input')

%% Fourier series representation of signal (Phase Spectrum)
phase = angle(c_in);
phase(1:floor(length(phase)/2)) = -phase(1:floor(length(phase)/2));

nexttile;
Hp1=stem(f,phase);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
axis([-8*f0 8*f0 -pi pi])
xlabel('Frequency (Hz)'); ylabel('Phase');
title('Phase Spectrum of Input')

%% 
figure(2);
output_plots = tiledlayout(2,2);
title(output_plots, 'Output Plots', 'FontSize', 20);

%% Designing the 2nd order Butterworth filter parameters
fc = 11500; % set your cutoff frequency

Hf = 1 ./ (1 + 1.414*(1i*f/fc) + (1i*f/fc).^2); % transfer function of filter
c_out = c_in .* Hf; % Fourier coefficients of the filter output

%% Construct the output signal from the Cout Fourier coefficients
A = zeros(2*N+1,ceil(no_sample));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*tt);
end
output = sum(A);
nexttile([1 2]);
Hp1 = plot(tt,input,'r',tt,real(output),'b');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
axis([tt(1) tt(length(tt)) min([min(real(output)) min(input)]) max([max(real(output)) max(input)])]);
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Filter Input and Output - Time Domain')
set(Ha,'Fontsize',16)
legend('Input','Output')

%% Frequency response of output

% Plot input and output magnitude spectrum
nexttile;
stem(f,abs(c_in),'r','LineWidth',2);
hold on
stem(f,abs(c_out),'b','LineWidth',2);
hold off
axis([-8*f0 8*f0 0 max(abs(c_in))])
Ha = gca;
set(Ha,'Fontsize',16)
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Magnitude Spectrum of Filter Output and Input')
Ha = gca;
set(Ha,'Fontsize',16)
legend('Input','Output')

%% Fourier series representation of signal (Phase Spectrum)
phase_output = angle(c_out);
phase_output(1:floor(length(phase_output)/2)) = phase_output(1:floor(length(phase_output)/2));

nexttile;
stem(f,phase,'r','LineWidth',2);
hold on
stem(f,phase_output,'b','LineWidth',2);
hold off
axis([-8*f0 8*f0 -pi pi])
Ha = gca;
set(Ha,'Fontsize',16)
xlabel('Frequency (Hz)'); ylabel('Phase');
title('Phase Spectrum of Filter Output and Input')
Ha = gca;
set(Ha,'Fontsize',16)
legend('Input','Output')

% %% Frequency response of output
% magFFT = abs(fft(output))/no_sample;
% fshift = (-no_sample/2:no_sample/2-1)*(fs/no_sample);
% 
% % plot amplitude spectrum
% figure(6)
% Hp1 = plot(fshift, fftshift(magFFT));
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% xlim([-8*f0 8*f0]);
% xlabel('Frequency (Hz)'); ylabel('Amplitude');
% title('Magnitude Spectrum of Output');