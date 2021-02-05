clc
clear all

%% Square wave generator (input)
f0 = 10000; % fundamental freq of input square wave
T0 = 1/f0;  % period
tstep = 0.001*T0;
tt = -3*T0:tstep:3*T0;
fs = 1/tstep;   % sampling frequency

input = square(tt*2*pi*f0,50); % input square wave
N = length(input);

% plot square wave
figure(1)
Hp1 = plot(tt,input);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
xlim([tt(1) tt(length(tt))]);
title('Input - Time Domain')

%% Frequency response of input
magFFT = abs(fft(input))/N;
fshift = (-N/2:N/2-1)*(fs/N);

% plot amplitude spectrum
figure(2)
Hp1 = plot(fshift, fftshift(magFFT));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
xlim([-8*f0 8*f0]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Magnitude Spectrum of Input');

% plot phase spectrum
phaseFFT = angle(fft(input));
figure(3)
Hp1 = plot(fshift, fftshift(phaseFFT));
xlim([-8*f0 8*f0]);

%% Designing the 2nd order Butterworth filter parameters
n = 2; % order of butterworth filter
fc = 11500; % set your cutoff frequency
fund_freq = 10e3; % fundamental frequency
third_harm_freq = 30e3; % third harmonic frequency

% a for poles (or factors from textbook)
switch n
    case 1
        a = [1 1]; % first order
    case 2
        a = [1 1.414 1]; % second order
    case 3
        a = conv([1 1], [1 1 1]); % third order
    case 4
        a = conv([1 0.765 1], [1 1.848 1]); % fourth order
end

% b for zeros
b = 1;

%% calculate normalized frequency response
% range = ceil(third_harm_freq/fc) + 1;
%w = linspace(-range,range,length(f));
%h = freqs(b,a,f);
%mag = 20*log10(abs(h)); %% convert magnitude to dB
%phase = angle(h);
%phasedeg = phase*180/pi;
Hf = 1 ./ (1 + 1.414*(1i*fshift/fc) + (1i*fshift/fc).^2);

output = ifft(fftshift(fft(input)) .* Hf);
figure(4)
Hp1 = plot(tt, abs(output),'b',tt,input,'r');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('filter input and output-time domain')
set(Ha,'Fontsize',16)
legend('output','input')
xlim([tt(1) tt(length(tt))])

%% Frequency response of output
magFFT_output = abs(fft(output))/N;

% plot amplitude spectrum
figure(5)
Hp1 = plot(fshift, fftshift(magFFT_output));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
xlim([-8*f0 8*f0]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Magnitude Spectrum of Output');

