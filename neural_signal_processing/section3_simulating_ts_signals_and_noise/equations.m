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
% optional prettification of the plot
%set(gca,'xlim',[-1.1 1.1],'ylim',[-1.1 1.1]); 
xlabel('Time (ms)')
title('My first sine wave plot! Mom will be so proud!')

%% the sum of sine waves can appear like a complicated time series

% create a few sine waves and sum

% define a sampling rate
srate = 1000;

% list some frequencies
frex = [ 3   10   5   15   35 ];

% list some random amplitudes... make sure there are the same number of
% amplitudes as there are frequencies!
amplit = [ 5   15   10   5   7 ];

% phases... list some random numbers between -pi and pi
phases = [  pi/7  pi/8  pi  pi/2  -pi/4 ];

% define time...
time = -1:1/srate:1;


% now we loop through frequencies and create sine waves
sine_waves = zeros(length(frex),length(time));
for fi=1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
end

% now plot the result
figure(2), clf
plot(time,sum(sine_waves))
title('sum of sine waves')
xlabel('Time (s)'), ylabel('Amplitude (arb. units)')


% now plot each wave separately
figure(3), clf
for fi=1:length(frex)
    subplot(length(frex),1,fi)
    plot(time, sine_waves(fi, :))
    axis([ time([1 end]) -max(amplit) max(amplit) ])
end