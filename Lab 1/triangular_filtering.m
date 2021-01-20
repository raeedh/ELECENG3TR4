%triangular wave filtering
clear all
hold off
j=sqrt(-1);
f0=500;     %fundamental freq of input triangular wave
T0 = 1/f0;  %period 
tstep = 0.005*T0;
no_sample = 3*T0/tstep + 1; %no. of samples  within  3*T0
no_sample1 = T0/tstep + 1; %no. of samples  within  T0
%tt = -0.5*T0:tstep:0.5*T0;
tt = -1.5*T0:tstep:1.5*T0;

tt1 = -0.5*T0:tstep:0.5*T0; % time vector for the period -0.5T0 to 0.5T0
gp1 = tt1/T0; %input - triangular wave in the period -0.5T0 to 0.5T0
gp_in = [gp1 gp1(2:no_sample1-1) gp1]; %3 cycles of the triangular wave
figure(1)
Hp1 = plot(tt,gp_in)
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('input - time domain')
pause

      
      
K=1/(2*pi);
N=100; %no. of harmonics
nvec = -N:N;
for n = nvec
    m = n+N+1;
    absn = abs(n);
    if(n > 0)
      c_in(m) = j*K*((-1)^n)/n;
    elseif (n == 0)
      c_in(m) = 0.0;
    else 
      c_in(m) = -j*K*((-1)^absn)/absn;
    end
end
f = nvec*f0; %frequency vector
figure(2)
Hp1=plot(f,abs(c_in))
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of input')
pause
figure(3)
Hp1=plot(f,angle(c_in))
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
axis([-0.1e4 0.1e4 -pi pi])
title('phase spectrum of input')
pause
R=3.3e3;
C=0.1e-6;
fc=1/(2*pi*R*C)     %cutoff freq of filter
%fc = 5000;
pause

Hf = 1 ./(1+j*f/fc) ;%filter transfer function
c_out = c_in .* Hf; %Fourier coefficients of the filter output
figure(4)
Hp1=plot(f,abs(c_out),'b',f,abs(c_in),'r')
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of filter output and input')
% pause(2)
% hold on
% Hp1=plot(f,abs(c_in),'r')
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
legend('output','input')
pause
% hold off
% Hp1=plot(f,angle(c_out))
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% title('phase spectrum of output')
% axis([-0.1e4 0.1e4 -pi pi])
% pause
% hold on
% Hp1=plot(f,angle(c_in),'r')
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% pause
% hold off

%c_out = c_in;

A = zeros(2*N+1,ceil(no_sample));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(j*2*pi*n*f0*tt);
end
gp_out = sum(A);
figure(5)
Hp1 = plot(tt,real(gp_out),'b',tt,gp_in,'r');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('filter input and output-time domain')
% pause
% hold on
% Hp1 = plot(tt,gp_in,'r')
% set(Hp1,'LineWidth',2)
% Ha = gca;
set(Ha,'Fontsize',16)
legend('output','input')