%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Simulating EEG data
%      VIDEO: Project 2: dipole-level EEG data
% Instructor: sincxpress.com
%
%%

%% 

% mat file containing EEG, leadfield and channel locations
load emptyEEG

% select dipole location (more-or-less random)
diploc = 109;

% plot brain dipoles
figure(1), clf, subplot(121)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'rs','markerfacecolor','k','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')


% Each dipole can be projected onto the scalp using the forward model. 
% The code below shows this projection from one dipole.
subplot(122)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole projection')

%% add signal to one dipole and project to scalp

% reduce data size a bit
EEG.pnts  = 2000;
EEG.times = (0:EEG.pnts-1)/EEG.srate;

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);

% add signal to one dipole
dipole_data(diploc,:) = sin(2*pi*10*EEG.times);

% now project dipole data to scalp electrodes
EEG.data = squeeze(lf.Gain(:,1,:))*dipole_data;

% plot the data
plot_simEEG(EEG,31,2);

%% now for the projects!

%%%% IMPORTANT! Check internal consistency with existing EEG structure!!!

EEG


%% 1) pure sine wave with amplitude explorations

EEG.trials = 40;
EEG.data = zeros(64, 2004);
EEG.pnts = 2004;
EEG.times = (0:EEG.pnts-1)/EEG.srate;

% dipole amplitude magnitude
ampl = 1;
freq = 10;

% initialize all dipole data
sinewave = ampl * sin(2 * pi * freq * EEG.times);

% compute one trial
EEG.data(1, :) = sinewave;
lf.Gain(1, 1, :) = sinewave;

diploc = 1;

% plot brain dipoles
figure(1), clf, subplot(121)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'rs','markerfacecolor','k','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')

% Each dipole can be projected onto the scalp using the forward model. 
% The code below shows this projection from one dipole.
subplot(122)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole projection')



% repeat that for N trials
for trial = 1:EEG.trials
    EEG.data(trial, :) = ampl * sin(2 * pi * freq * EEG.times);
end

% plot the data
plot_simEEG(EEG,1,3);
plot_simEEG(EEG, 15, 4)

%%% Question: What is the smallest amplitude of dipole signal that still
%             elicits a scalp-level response?

%% 2) sine wave with noise

%%% Question: Given amplitude=1 of dipole signal, what standard deviation of noise
%             at all other dipoles overpowers the signal (qualitatively)?

% noise standard deviation
noiseSD = 5;
freq = 10;

for trial = 1:EEG.trials
    dipole_data = noiseSD * randn(size(lf.Gain, 3), EEG.pnts);
    dipole_data(diploc, :) = sin(2*pi*freq*EEG.times);

    EEG.data(:, :, trial) = squeeze(lf.Gain(:, 1, :)) * dipole_data;
end

plot_simEEG(EEG, 31, 2)





% % initialize all dipole data
% sinewave = ampl * sin(2 * pi * freq * EEG.times);
% 
% % compute one trial
% EEG.data(1, :) = sinewave +  2 * pi * randn(size(EEG.times));
% lf.Gain(1, 1, :) = sinewave + 2 * pi * randn(size(EEG.times));
% 
% diploc = 1;
% 
% % plot brain dipoles
% figure(1), clf, subplot(121)
% plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
% hold on
% plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'rs','markerfacecolor','k','markersize',10)
% rotate3d on, axis square
% title('Brain dipole locations')
% 
% % Each dipole can be projected onto the scalp using the forward model. 
% % The code below shows this projection from one dipole.
% subplot(122)
% topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
% set(gca,'clim',[-1 1]*40)
% title('Signal dipole projection')
% colormap jet
% colorbar
% 
% % plot the data
% plot_simEEG(EEG,1,2);
% colormap jet
% colorbar


%% 3) Non-oscillatory transient in one dipole, noise in all other dipoles

% simulation parameters
ptime = 1;  % peak time 
width  = .12;
ampl = 70;

gaus = ampl * exp(-(EEG.times - ptime).^2 / (2*width^2));

for trial = 1: EEG.trials
    dipole_data = randn(size(lf.Gain, 3), EEG.pnts);
    dipole_data(diploc, :) = gaus;

    EEG.data(:, :, trial) = squeeze(lf.Gain(:, 1, :)) * dipole_data;
end

plot_simEEG(EEG, 31, 2)

% % Gaussian
% gwin = ampl * exp( -(4*log(2) * (EEG.times-ptime).^2) / fwhm^2);
% 
% % empirical FWHM
% gwinN   = gwin./max(gwin);
% midp    = dsearchn(EEG.times', ptime);
% pst5    = midp-1+dsearchn(gwinN(midp:end)',.5);
% pre5    = dsearchn(gwinN(1:midp)',.5);
% empfwhm = EEG.times(pst5) - EEG.times(pre5);
% 
% 
% figure(3), clf, hold on
% plot(EEG.times,gwin,'k','linew',2)
% plot(EEG.times([pre5 pst5]),gwin([pre5 pst5]),'ro--','markerfacecolor','k')
% plot(EEG.times([pre5 pre5]),[0 gwin(pre5)],'r:')
% plot(EEG.times([pst5 pst5]),[0 gwin(pst5)],'r:')
% title([ 'Requested FWHM: ' num2str(fwhm) 's, empirical FWHM: ' num2str(empfwhm) 's' ])
% xlabel('Time (s)'), ylabel('Amplitude')
% 
% 
% 
% EEG.data(1, :) = gwin;
% lf.Gain(1, 1, :) = gwin;
% 
% % noise in other dipoles
% hz = linspace(0, EEG.srate/2, floor(length(EEG.times) / 2) - 1);
% 
% for trial = 2:EEG.trials
%     % generate 1/f amplitude
%     ed = 50; % exponential decay -> defines how the amplitude function goes down
%     as = rand(1, floor(EEG.pnts / 2) - 1) .* exp(-(1:floor(EEG.pnts/2) - 1) / ed);
%     as = [as(1) as 0 0 as(:, end:-1:1)];
%     
%     % Fourier coefficient
%     fc = as .* exp(1i * 2 * pi * rand(size(as)));
%     
%     % inverse Fourier transform to create noise
%     noise = real(ifft(fc)) * EEG.pnts;
% 
%     EEG.data(trial, :) = noise(1, 1:end-1);
%     lf.Gain(trial, 1, :) = noise(1, 1:end-1);
% end
% 
% figure(4)
% topoplotIndie(-lf.Gain(:,1,1), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
% set(gca,'clim',[-1 1]*40)
% title('Signal dipole projection')
% 
% 
% figure(5)
% topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
% set(gca,'clim',[-1 1]*40)
% title('Signal dipole projection')
% 
% plot_simEEG(EEG,31,2);

%% 4) Non-stationary oscillation in one dipole, transient oscillation in another dipole, noise in all dipoles

%%% first pick two dipoles
dipole1 = 15;
dipole2 = 37;

% noise in all dipoles
for t = 1:length(EEG.trials)
    for d = 1:length(EEG.nbchan)
        EEG.data(d, :, t) = 3 * randn(size(EEG.times)) + 5;
    end
end

% non-stationary oscillation
for i = 1:length(EEG.trials)
    non_stationary = sin(2*pi*10*EEG.times + 5 * pi/6 * randn);
    EEG.data(dipole1, :, i) = non_stationary;
end 

% transient oscillations
peaktime = 1; % seconds
width = .12;

for t = 1:length(EEG.trials)
    gaus = exp(-(EEG.times - peaktime).^2 / (2*width^2));
    sinewave = sin(2*pi*10 * EEG.times ) .* gaus;
    EEG.data(dipole2, :, t) = sinewave;
end


%%% then do the simulation


% plot the data
plot_simEEG(EEG,dipole1,3);
plot_simEEG(EEG,dipole2,2);


%% 
