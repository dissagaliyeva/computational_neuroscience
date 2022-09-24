% y_t = sum(x_i, g_i)_{i=t-k}^(t+k)

%%
srate = 3000; % Hz
time = linspace(1, srate + 1, srate + 1);
n = length(time);
p = 15;

noiseamp = 5;
amp = interp1(rand(p, 1) * 30, linspace(1, p, n));
noise = noiseamp * randn(size(time));
signal = amp + noise;

%%
plot(time, signal);
hold on
plot(time, amp)

%% create Gaussian kernel

% full-width half-maximum: the key Gaussian parameter
fwhm = 25; % in ms

% normalied time vector in ms
k = 100;
gtime = 1000 * (-k : k) / srate;

% create Gaussian window
gauswin = exp((-4 * log(2) * gtime.^2) / fwhm^2);

% compute empirical FWHM
prePeakHalf = k + dsearchn(gauswin(k + 1 : end)', .5);
pstPeakHalf = dsearchn(gauswin(1:k)', .5);

empFWHM = gtime(prePeakHalf) - gtime(pstPeakHalf);

% show Gaussian
figure(1), clf, hold on
plot(gtime, gauswin, 'ko-', 'markerfacecolor', 'w', 'linew', 2);
plot(gtime([prePeakHalf, pstPeakHalf]), gauswin([prePeakHalf, pstPeakHalf]), 'm', 'linew', 3)

% normalize Gaussian to unit energy
gauswin = gauswin / sum(gauswin);
title(['Gaussian kernel with requetes FWHM ' num2str(fwhm) ' ms '])
xlabel('Time (ms)'), ylabel('gain')

%% implement the filter

% intialize filtered signal vector
filtsigG = signal;

% implement the running mean filter
for i = k + 1 : n - k - 1
    % for each point is the weighted average of k surrounding points
    filtsigG(i) = sum(signal(i - k : i + k).*gauswin);
end

plot(2), clf, hold on
plot(time, signal, 'r')
plot(time, filtsigG, 'k', 'linew', 3)
xlabel('Time (s)'), ylabel('amp. (a.u.)')
legend({'Original signal'; 'Gaussian-filtered'});
title('Gaussian smoothing filter');

