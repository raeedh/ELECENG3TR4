clc
clear all

%% Square wave generator (input)
f0 = 10000; % fundamental freq of input square wave
T0 = 1/f0;  % period
tstep = 0.005*T0;
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
inputFFT = abs(fft(input))/N;
fshift = (-N/2:N/2-1)*(fs/N);

% plot amplitude spectrum
figure(2)
Hp1 = plot(fshift, fftshift(inputFFT));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
xlim([-8*f0 8*f0]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Magnitude Spectrum of Input');


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
%Hf = 1 ./ (1 + 1.414*(1i*f/fc) + (1i*f/fc).^2); 

%c_out = c_in .* Hf ; %Fourier coefficients of the filter output

% figure(4)
% stem(f,abs(c_in),'r','LineWidth',2);
% hold on
% stem(f,abs(c_out),'b','LineWidth',2);
% hold off
% axis([-8*f0 8*f0 0 max(abs(c_in))])
% Ha = gca;
% set(Ha,'Fontsize',16)
% title('magnitude spectrum of filter output and input')
% Ha = gca;
% set(Ha,'Fontsize',16)
% legend('input','output')
% pause

%% Construct the output signal from the Cout Fourier coefficients

% A = zeros(2*N+1,ceil(no_sample));
% for n = nvec
%     m=n+N+1;
%     A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*tt);
% end
% gp_out = sum(A);
% figure(5)
% Hp1 = plot(tt,real(gp_out),'b',tt,gp1,'r');
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% title('filter input and output-time domain')
% set(Ha,'Fontsize',16)
% legend('output','input')
