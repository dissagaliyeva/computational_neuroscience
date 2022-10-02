%% Sine waves

% define variables
freq = 2; % freq in Hz
srate = 1000; % (1 second)
time = -1:1/srate:1;
ampl = 2;
phase = pi/3;

% create sinewave
sinewave = ampl * sin(2*pi*freq*time + phase);

plot(time, sinewave)