%% filter parameters
n = 1; % order of butterworth filter

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