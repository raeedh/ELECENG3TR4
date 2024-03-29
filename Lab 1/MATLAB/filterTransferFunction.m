clear;
clc;

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
cont_trans_func = tf(b, a)
dis_trans_func = c2d(cont_trans_func, 1e-5)

range = ceil(third_harm_freq/fc) + 1;
w = linspace(-range,range,500);
[mag, phase, wout] = bode(dis_trans_func, w);

magdb = 20 * log10(mag(:));
phase = phase(:);

%% plot frequency response
% normalized frequency response
figure(10)
subplot(2,1,1)
plt1 = plot(w,magdb);
xlim([w(1) w(length(w))]);
grid on
title('Normalized Frequency Response (Magnitude)');
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
xline(1,'k',{'Cutoff frequency'}); % 3dB frequency (at cutoff frequency)


subplot(2,1,2)
plt2 = plot(fc*w,magdb);
xlim([fc*w(1) fc*w(length(w))]);
plt2ax = ancestor(plt2(1), 'axes'); plt2ax.XAxis.Exponent = 0;
grid on
title('Frequency Response (Magnitude)');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xline(fc,'b',{'Cutoff frequency'}); % 3dB frequency (at cutoff frequency)
xline(fund_freq,'r',{'Fundamental frequency'}); % first harmonic
xline(third_harm_freq,'m',{'3rd harmonic'}); % third harmonic

set(gcf, 'WindowState', 'maximized');

% phase response
%plot(w,phasedeg)
%grid on
%xlabel('Frequency (rad/s)')
%ylabel('Phase (degrees)')

%% determine gains at desired frequencies
vq = interp1(fc*w, magdb, [fc fund_freq third_harm_freq]);

message = ['Attenuation at cutoff frequency: ', num2str(vq(1)), ' dB'];
disp(message);

message = ['Attenuation at fundamental frequency: ', num2str(vq(2)), ' dB'];
disp(message);

message = ['Attenuation at third harmonic: ', num2str(vq(3)), ' dB'];
disp(message);

switch (abs(vq(2)) < 2) 
    case 1
        temp = 'Yes';
    case 0
        temp = 'No';
end
message = ['The loss of the fundamental frequency due to filtering is less than 2 dB: ', temp];
disp(message);

message = ['Output of harmonics is at least ', num2str(-(vq(3)-vq(2))), ' dB below fundamental frequency'];
disp(message);

switch (abs(vq(3)-vq(2)) > 13.46) 
    case 1
        temp = 'Yes';
    case 0
        temp = 'No';
end
message = ['This meets the required attentuation from filtering (at least 13.46 dB): ', temp];
disp(message);
