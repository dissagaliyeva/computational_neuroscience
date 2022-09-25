%% pre-reqs
% before jumping to Fourier transform, we need to know about
% sine waves, dot product, and complex numbers

% sine wave example
srate = 1000;  % sampling rate in Hz (1s)
time = 0:1/srate:5; % units of seconds (from 0 to 5 with step size 1/1000)
f = 5; % frequency, units in Hz
a = 1; % amplitude, arbitrary unit
theta = pi / 2; % in radians
sinewave = a * sin(2 * pi * f * time + theta);

plot(time, sinewave)

%% complex numbers + sine waves

csw = exp(1i * 2 * pi * f * time);
plot3(time, real(csw), imag(csw))


%% Slow Fourier Transform

N = length(time);
signalX = zeros(N);
fTime = (0 : N- 1) / N;

for fi=1:N
    fSine = exp(-1i * 2 * pi * (fi - 1).*fTime);
    signalX(fi) = sum(fSine.*signalX);
end

signalX = signalX / N;

%% getting lowest and highest frequencies 

low = 0;
high = srate / 2;  % nyquist freq

% next, we need to know how many frequency steps in the range [0, nyquist
% freq) -> freq. resolution

hz = linspace(0, high, N/2 + 1);
% hz = (0:1/srate:N/2+1) * high;   % alternative





