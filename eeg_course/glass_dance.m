% load the file
load data\glassDance.mat

% play the music
% soundsc(glassclip, srate)

% some variables for convenience
pnts = length(glassclip);
timevec = (0:pnts - 1) / srate;

% draw the time-domain signals
figure(1), clf
subplot(511)
plot(timevec, glassclip)
xlabel('Time (s)')

% static power spectrum and pick a frequency range

% inspect power spectrum
hz = linspace(0, srate / 2, floor(length(glassclip) / 2) + 1);
powr = abs(fft(glassclip(:, 1)) / pnts);

subplot(512), cla 
plot(hz, powr(1:length(hz)))
set(gca, 'xlim', [100 2000], 'ylim', [0 max(powr)])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

% pick freqs to filter
frange = [300 460];

% design an FIR1 filter
fkern = fir1(2001, frange / (srate / 2), 'bandpass');

% apply the filter to the signal
filtglass(:, 1) = filtfilt(fkern, 1, glassclip(:, 1));
filtglass(:, 2) = filtfilt(fkern, 1, glassclip(:, 2));

% plot filtered signal power spectrum
hold on
powr = abs(fft(filtglass(:, 1)) / pnts);
plot(hz, powr(1:length(hz)), 'r')

% plot time-frequency response (spectrogram)
subplot(5, 1, 3:5)
spectrogram(glassclip(:, 1), hann(round(srate / 10)), [], [], srate, 'yaxis')
hold on 
plot(timevec([1 1; end end]), frange([1 2; 1 2]) / 1000, 'k:', 'linew', 2)

set(gca, 'ylim', [0, 2])

% play the sound
soundsc(filtglass, srate)

