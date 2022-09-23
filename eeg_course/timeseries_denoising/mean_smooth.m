% Running-mean time series filter
% y_t = (2k + 1)**-1 * sum(x_i)_{i=t-k}^{t+k}

%%

% create signal
srate = 1000; % Hz
time = 0:1 / srate:3;
n = length(time);
p = 15; % poles for random interpolation

% noise level, measured in std
noiseamp = 5;

% amplitude modulator and noise level
amp = interp1(rand(p, 1) * 30, linspace(1, p, n));
noise = noiseamp * randn(size(time));
signal = amp + noise;

%%

plot(time, signal, 'r')
hold on 
plot(time, amp, 'b')

%%
% initialize filtered signal vector
% setting filtsig=signal will remove the edge points problem and set it to
% the starting/ending points instead
filtsig = zeros(size(signal));

% implement running mean filter
% making k bigger will give smoother lines
k = 20;

for i=k+1:n-k-1
    % each point is the average of k surrounding points
    k
    filtsig(i) = mean(signal(i - k: i + k));
end

% compute window size in ms
windowsize = 1000 * (k * 2 + 1) / srate;

% plot the noisy and filtered signals
figure(1), clf, hold on
plot(time, signal, time, filtsig, 'linew', 2)

% draw patch to indicate the window size
tidx = dsearchn(time', 1);
ylim = get(gca, 'ylim');
patch(time([tidx - k tidx - k tidx + k tidx + k]), ylim([1 2 2 1]), 'k', 'facealpha', .2)
plot(time([tidx tidx]), ylim, 'k--')


