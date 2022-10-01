%%
%   COURSE: Neural signal processing and analysis: Zero to hero
%  SESSION: Problem set: Simulating and visualizing data
%  TEACHER: Mike X Cohen, sincxpress.com
%

%%

%%% INSTRUCTIONS:
% The goal of this assignment is to simulate time series data
% that can be used to test time-series analysis methods.
% For each section below: 
%   1) Complete the MATLAB code
%   2) Put the data into the EEG structure
%      - Make sure all relevant fields are accurate (EEG.data, EEG.pnts, EEG.trials, EEG.srate, EEG.nbchan, EEG.times)
%   3) Plot the data using the function plot_simEEG()

% NOTE 1: Obviously, you need to fill in missing code.
% NOTE 2: Be careful, because sometimes there is incorrect code that doesn't produce coding errors.
%         Remember: Visualize, and visualize often.

%% 1) white noise

% The goal of this exercise is to gain basic familiarity with data simulations.
% You will create a dataset of white noise and plot it.

% specify EEG parameters
EEG.srate  = 500; % sampling rate in Hz
EEG.pnts   = EEG.srate * 5; % doing so will automatically convert it to seconds (here 2s)
EEG.trials = 20;
EEG.nbchan = 20;

% time vector
EEG.times = (0:EEG.pnts-1)/EEG.srate;

% create data as white noise
EEG.data = randn(EEG.nbchan, EEG.pnts, EEG.trials);

% the function below takes at least one argument (EEG),
% and optionally a second argument (channel number),
% and optionally a third argument (figure number)
plot_simEEG(EEG,15,1)
colormap jet
colorbar

%%% Question: What is the effect of noise amplitude on the 
%             resulting graphs?

EEG.data = 10 * randn(EEG.nbchan, EEG.pnts, EEG.trials);
plot_simEEG(EEG,15,1)
colormap jet
colorbar

%%% Question: Do the results change if you use normally distributed
%             vs. uniformly distributed noise?

EEG.data = rand(EEG.nbchan, EEG.pnts, EEG.trials);
plot_simEEG(EEG,2,1)
colormap jet
colorbar

set(gca, 'ylim', [0 .001])  % run this only by F9 after the plots are shown


% this uniform distribution has a zero-component that drags the whole
% figure, therefore, to remove it run the following:

EEG.data = rand(EEG.nbchan, EEG.pnts, EEG.trials) - 0.5;
plot_simEEG(EEG,2,1)
colormap jet
colorbar

%%% Question: Are the results different for different channels? Why or why not?


%% 2) pink noise

% The goal of this exercise is to extend the previous exercise to "pink" noise.
% You should create the noise separately on each trial


% feel free to change some parameters compared to above...
EEG.nbchan = 64;

% the key parameter of pink noise is the exponential decay (ed)
ed = 10; % try different values!

% initialize EEG data as a zeros matrix
EEG.data = zeros(EEG.trials, EEG.pnts, EEG.nbchan);


for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        %%% Note about the code below: This involves creating the signal first in the frequency domain, 
        %   then transforming it into the time domain. Don't worry if you don't understand the details
        %   (you'll learn them tomorrow!); try to plot each step to get a basic idea.
        
        % generate one-sided 1/f amplitude spectrum
        as = rand(1,EEG.pnts) .* exp(-(0:EEG.pnts-1)/ed);
        
        % Fourier coefficients as amplitudes times random phases
        fc = as .* exp(1i*2*pi*rand(size(as)));
        
        % inverse Fourier transform to create the noise
        EEG.data(chani,:,triali) = real(ifft(fc));
    end
end

plot_simEEG(EEG,5,1)
colormap jet
colorbar

%%% Question: Which looks more like real EEG data: white or pink noise?
%             Why do you think this is?
% 
%%% Question: Which values of variable 'ed' make the data look most like real EEG data?


%% 3) Ongoing stationary

% The goal here is to create a dataset with ongoing sinewaves.
% There should be multiple sine waves simultaneously in each channel/trial.

% list of frequencies and corresponding amplitudes
frex = [ 3 5 16 ]; % in Hz
amps = [ 3 4 5  ]; % in arbitrary units

EEG.trials = 100;

% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        %%% Note that here the signal is created in the time domain, unlike in the previous example.
        %   Some signals are easier to create in the time domain; others in the frequency domain.
        
        % create a multicomponent sine wave
        sinewave = zeros(1,EEG.pnts);
        for si=1:length(frex)
            sinewave = sinewave + amps(si) * sin(2*pi*frex(si) * EEG.times);
        end
        
        % data as a sine wave plus noise
        EEG.data(chani, :, triali) = sinewave + randn(size(sinewave));
    end
end

% and plot
plot_simEEG(EEG, 1, 1)
colormap jet
colorbar

% in this case we see that in the time-frequency plot there are bursts at
% fruqncy 3-5; that is the articats in the tf domain. The real results
% should be a static line as in 13Hz instead of the bursts

%%

%%% Question: What can you change in the code above to make the EEG
%             activity non-phase-locked over trials?

% list of frequencies and corresponding amplitudes
frex = [ 3 5 16 ]; % in Hz
amps = [ 3 4 5  ]; % in arbitrary units

EEG.trials = 100;   

% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        %%% Note that here the signal is created in the time domain, unlike in the previous example.
        %   Some signals are easier to create in the time domain; others in the frequency domain.
        
        % create a multicomponent sine wave
        sinewave = zeros(1,EEG.pnts);

        sinewave = sinewave + amps(2) * sin(2*pi*(frex(2)) * EEG.times * rand * 2 * pi);
       
        % data as a sine wave plus noise
        EEG.data(chani, :, triali) = sinewave + randn(size(sinewave));
    end
end

% and plot
plot_simEEG(EEG, 1, 2)

% randi will work, but that function generates integers, so the range of 
% phase values is severely restricted. rand*2*pi will sample all possible phase values.


%%
%%% Question: Which of the plots look different for phase-locked vs. non-phase-locked?
%             (Hint: plot them in different figures to facilitate comparison.)
%             Are you surprised about the differences?
%             
%%% Question: Are all frequencies equally well represented in the 'static' and 'dynamic' power spectra?
%             Can you change the parameters to make the spectral peaks more or less visible in the two plots?
% 


%% 4) ongoing nonstationary

% Here you want to create narrowband non-stationary data. 
% This is starting to be more "realistic" (in a signal-characteristic sense) for EEG data.

% signal parameters in Hz
peakfreq = 14;
fwhm     =  1;

% frequencies
hz = linspace(0,EEG.srate,EEG.pnts);

%%% create frequency-domain Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-peakfreq;          % shifted frequencies
fg = exp(-.5*(x/s).^2);    % gaussian


for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        %%% As with previous simulations, don't worry if you don't understand the mechanisms;
        %   that will be clear tomorrow. Instead, you can plot each step to try to build intuition.
        
        % Fourier coefficients of random spectrum
        fc = rand(1,EEG.pnts) .* exp(1i*2*pi*rand(1,EEG.pnts));
        
        % taper Fourier coefficients by the Gaussian
        fc = fc .* fg;
        
        % go back to time domain to get EEG data
        EEG.data(chani,:,triali) = real( ifft(fc) );
    end
end

% plot magnitude coefficient of the Fourier with respect to amplitude
% spectrum 
% plot(hz, abs(fc.*fg))

plot(real(ifft(fc)))

% and plot
plot_simEEG(EEG, 1, 1)
plot_simEEG(EEG, 1, 2)
%             
%%% Question: What is the effect of FWHM on the results? Is larger or smaller more realistic?
% the resulting fwhm looks best and most natural when it's equal to 5.
% There should be a distinctive tf plot with bold line; 
% 
%%% Question: Can you modify the code to have narrowband activity at two different frequency ranges?
% 

% signal parameters in Hz
peakfreq = 32;
peakfreq2 = 15;
fwhm     =  5;

% frequencies
hz = linspace(0,EEG.srate,EEG.pnts);

%%% create frequency-domain Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-peakfreq;          % shifted frequencies
fg1 = exp(-.5*(x/s).^2);    % gaussian

%%% create frequency-domain Gaussian
s2  = fwhm*(2*pi-1)/(4*pi); % normalized width
x2  = hz-peakfreq2;          % shifted frequencies
fg2 = exp(-.5*(x2/s2).^2);    % gaussian

for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        %%% As with previous simulations, don't worry if you don't understand the mechanisms;
        %   that will be clear tomorrow. Instead, you can plot each step to try to build intuition.
        
        % Fourier coefficients of random spectrum
        fc = rand(1,EEG.pnts) .* exp(1i*2*pi*rand(1,EEG.pnts));
        
        % taper Fourier coefficients by the Gaussian
        fc1 = fc .* fg1;

        % taper Fourier coefficients by the Gaussian
        fc2 = fc .* fg2;
        
        % go back to time domain to get EEG data
        EEG.data(chani,:,triali) = real( ifft(fc1) ) + real(ifft(fc2));
    end
end

plot_simEEG(EEG)

%% 5) transients #1: Gaussian

% All the exercises above were for ongoing signals. Now we move to transients.
% Start with a Gaussian.

% gaussian parameters (in seconds)
peaktime = 1;
width = .1;

% re-initialize EEG data
EEG.data = zeros(EEG.nbchan,EEG.pnts,EEG.trials);

for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        % generate time-domain gaussian: normal
        % gaus = exp(-(EEG.times - peaktime).^2 / (width^2*2));
        
        % with jitter
        trialpeak = peaktime + randn * 100;
        gaus = exp(-(EEG.times - trialpeak).^2 / (width^2*2));

        % data are the Gaussian
        EEG.data(chani,:,triali) = gaus;
    end
end

% and plot
plot_simEEG(EEG,1,2);

%%% Questions: What happens if you add random jitter to the peaktime on each trial? 
%              How much jitter until the ERP is nearly gone?


%% 6) transients #2: oscillations w/ Gaussian

% Finally, we get to the most useful simulations for time-frequency analysis:
%   time-limited narrow-band activity. This is done by multiplying a Gaussian with a sine wave.

% sine wave frequency
sfreq = 8;

% gaussian parameters (in seconds)
peaktime = 1;
width = .2;


% re-initialize EEG data
EEG.data = zeros(EEG.nbchan,EEG.pnts,EEG.trials);

for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        % generate time-domain gaussian
        trialpeak = peaktime + randn/5;
        gaus = exp( -(EEG.times-trialpeak).^2 / (2*width^2) );
        
        
        % generate sine wave with same phase
        sw = cos(2*pi*sfreq*EEG.times);
        
        
        % data are sine wave times Gaussian
        EEG.data(chani,:,triali) = gaus;
    end
end

% and plot
plot_simEEG(EEG,1,1);

%             
%%% Question: Do the results look realistic? What can you change to make it look even more EEG-like?
% 
% 
%             
%%% Question: How can you modify the code to make the transient non-phase-locked?
%             Which of the three data plots are most affected by phase-locked vs. non-phase-locked?



%%


%% More exercises for more fun!


%% 7) Add pink noise to #5


%% 8) Combine #4 and #6 to make a transient (Gaussian-windowed) non-stationary signal


%% done.
