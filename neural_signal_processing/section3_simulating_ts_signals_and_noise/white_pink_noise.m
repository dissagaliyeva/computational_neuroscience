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


%% Pink noise (fractal)

% simulation details 
srate = 500;
time = -1 : 1/srate : 5;
pnts = length(time);
hz = linspace(0, srate/2, floor(length(time) / 2) - 1);

% generate 1/f amplitude
ed = 50; % exponential decay -> defines how the amplitude function goes down
as = rand(1, floor(pnts / 2) - 1) .* exp(-(1:floor(pnts/2) - 1) / ed);
as = [as(1) as 0 0 as(:, end:-1:1)];

% Fourier coefficient
fc = as .* exp(1i * 2 * pi * rand(size(as)));

% inverse Fourier transform to create noise
noise = real(ifft(fc)) * pnts;

figure(5), clf
subplot(211)
plot(time, noise, 'k')
set(gca, 'fontsize', 12)
title('Pink noise: Time domain')
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
