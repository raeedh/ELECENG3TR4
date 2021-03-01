clear
hold off
format long e
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
% Am=1;
fm = 1e3;
% mt = Am*cos(2*pi*fm*tt);
Tm = 0.0005;
mt = -2*sinc(tt/Tm);

%max of absolute of m(t)
% maxmt = Am;
maxmt = max(abs(mt));
%For 50% modulation
ka=0.5/maxmt;

%AM signal
st = (1+ka*mt).*ct;

figure(1)
Hp1 = plot(tt,ct);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('Carrier c(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Carrier : Time domain');
axis([-1e-3 1e-3 -1.1 1.1])
pause(1)
figure(2)
Hp1 = plot(tt,mt);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('message  m(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('message signal : Time domain');
axis([-2e-3 2e-3 min(mt) max(mt)])
pause(1)
figure(3)
Hp1 = plot(tt,st);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('modulated wave : Time domain');
axis([-2e-3 2e-3 min(st) max(st)])
pause(1)
Mf1 = fft(fftshift(mt));
Mf = fftshift(Mf1);
figure(4)
Hp1=plot(freq,abs(Mf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the message signal');
axis ([-5e3 5e3 0 max(abs(Mf))])
%pause(5)

Sf1 = fft(fftshift(st));
Sf = fftshift(Sf1);
figure(5)
Hp1=plot(freq,abs(Sf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the modulated wave');
axis ([-25e3 25e3 0 max(abs(Sf))])
%pause(5)

%Demodulation

%tolerance should be decided by the trial and error 
%so that the capacitor discharging instant (i.e when
% yt(n-1) exceeds st(n) by a small amount 'tole'
tole = 0.1;
%time constant RC
%This should be optimized to avoid envelope distortion 
% RC = 0.5*(1/fc + 1/fm);
% 2i
% RC = 1/fc;
% 2ii
% RC = 10*Tm;
% 2iii
RC = 1.2*Tm;
%Envelope detector
yt = zeros(1,N);
yt(1) = st(1);
n=1;
for t=tt
    if(n > 1)
     if(st(n) > yt(n-1))
         yt(n) = st(n);
     else
         if((yt(n-1)-st(n)) < tole)
           yt0 = yt(n-1);
           yt(n) = yt0;
           %time when C starts discharging
           tc = t;
         else
           yt(n) = yt0*exp(-(t-tc)/RC);
         end
     end
    end
    n=n+1;
end
yt(1)=yt(2);

figure(6)
Hp1 = plot(tt,yt);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('y(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('After the envelope detector');
axis([-2e-3 2e-3 0 max(yt)])
pause(1)
figure(7)
%DC removal and division by ka
yt1 = (yt - 1)/ka;
Hp1 = plot(tt,yt1,'r',tt,mt,'k');
legend('after DC removal','message signal')
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('y1(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('After the DC removal');
axis([-2e-3 2e-3 min(mt) max(mt)])
%pause
%Low pass filter to remove the ripple
%choose the cutoff frequency of the filter to be slightly higher than
%the highest freq of the message signal
% f0 = 1.1*fm;
% mt1 = rect_filt(yt1,freq,f0);
% figure(8)
% Hp1 = plot(tt,mt1);
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% Hx=xlabel('Time (sec) ');
% set(Hx,'FontWeight','bold','Fontsize',16)
% Hx=ylabel('m1(t)  (Volt)');
% set(Hx,'FontWeight','bold','Fontsize',16)
% title('After the low pass filter');
% axis([-2e-3 2e-3 min(mt1) max(mt1)])








     
