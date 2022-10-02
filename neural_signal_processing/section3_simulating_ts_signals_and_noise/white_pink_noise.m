%% Normally distributed noise 

% simulation details
srate = 100; % sampling rate in Hz
time = -1 : 1/srate : 2;
pnts = length(time); % 2 seconds 

% frequencies for the power spectrum
hz = linspace(0, srate / 2, floor(length(time) / 2) - 1);

% noise parameters
stretch = .5;
shift = 3;

% seed
% rng(42)

% generate random data
noise = stretch * randn(size(time)) + shift;

figure(4), clf
subplot(211)
plot(time, noise, 'k')
set(gca, 'fontsize', 12)
title('Normally distributed: Time domain')
xlabel('Time (s)'), ylabel('Amplitude')

subplot(223)
[y, x] = hist(noise, 100);
plot(x, y, 'k', 'linew', 2)
xlabel('Values'), ylabel('N per bin')
title('Signal histogram (distribution)')
set(gca, 'fontsize', 12, 'xlim', [min(x), max(x)])

subplot(224)
amp = abs(fft(noise) / pnts);
amp(2:end) = 2 * amp(2:end);
plot(hz, amp(1:length(hz)), 'k')
title('Frequency domain')
set(gca, 'fontsize', 12)