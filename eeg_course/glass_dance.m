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





