%% a for poles (or factors from textbook)
a0 = [1 0.765 1];
a1 = [1 1.848 1];
a = conv(a0, a1);

%% b for zeros
b = 1;

%% calculate normalized frequency response
w = linspace(-3,3);
h = freqs(b,a,w);
mag = 20*log10(abs(h)); %% convert magnitude to dB
%phase = angle(h);
%phasedeg = phase*180/pi;

%% plot frequency response
%subplot(2,1,1)
plot(w,mag)
grid on
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
xline(1); %% 3dB frequency (at cutoff frequency)

%subplot(2,1,2)
%plot(w,phasedeg)
%grid on
%xlabel('Frequency (rad/s)')
%ylabel('Phase (degrees)')