%% square wave generator
clc
clear all
hold off

f0= 10000;     %fundamental freq of input square wave, 0.1ms period
T0 = 1/f0;  %period 
tstep = 0.005*T0;
no_sample = 3*T0/tstep + 1; %no. of samples  within  3*T0
no_sample1 = int32(T0/tstep + 1); %no. of samples  within  T0
%tt = -0.5*T0:tstep:0.5*T0;
tt = -1.5*T0:tstep:1.5*T0;

tt1 = -0.5*T0:tstep:0.5*T0; % time vector for the period -0.5T0 to 0.5T0
gp1 = [zeros(1,floor(0.5*length(tt1))) ones(1,ceil(0.5*length(tt1)))]; %input - triangular wave in the period -0.5T0 to 0.5T0
gp_in = [gp1 gp1(2:no_sample1-1) gp1]; %3 cycles of the triangular wave
figure(1)
Hp1 = plot(tt,gp_in);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('input - time domain')
pause

%% Fourier series representation of signal (Amplitude Spectrum)
      
K=1/(2*pi);
N=100; %no. of harmonics
nvec = -N:N;
c_in = zeros(size(nvec));
for n = nvec
    m = n+N+1;
    c_in(m) = 1i*K*((-1)^n)/n;
    
    if (n == 0)
      c_in(m) = 0.0;
    end
end
f = nvec*f0; %frequency vector
figure(2)
Hp1=stem(f,abs(c_in));
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of input')
pause
