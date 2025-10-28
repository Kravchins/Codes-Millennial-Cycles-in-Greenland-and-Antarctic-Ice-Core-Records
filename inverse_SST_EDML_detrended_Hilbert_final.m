%% MATLAB Code for Millennial-Scale Spectral Analysis of Ice Core Data
%% Any use of this code must refer to
% Kravchinsky,V.A., Zhang,R., Borowiecki,R., Goguitchaichvili,A., Czarnecki,J., Czarnecki,A., Boers,N., Berger,A., & van der Baan,M. (2025)
% Millennial Cycles in Greenland and Antarctic Ice Core Records: Evidence of Astronomical Influence on Global Climate
% Journal of Geophysical Research: Atmospheres, 130(7), e2024JD042810. 
% https://doi.org/10.1029/2024JD042810
% Questions/suggestions to vadim@ualberta.ca

% Initialization and Data Import
close all; clear; clc;  % Clear workspace and close figures

% Load borehole data (Time vs. Parameter)
data = xlsread('EDML_Antarctica_example.xlsx');     %Load borehole data curve
age = data(:,1);      %Time
parameter = data(:,2);       %Parameter
% st1 = smooth(age,parameter,0.002,'loess'); %smooth data
st1 = parameter; % do NOT smooth data
ts1 = timeseries(st1,age,'Name','Reference'); % create the time series

% Define resampling interval and apply linear interpolation and
% normalization
% time=[0.02:0.02:150.28]'; % min time : interpolation interval : max time
time=[0.05:0.05:150.0]'; % min time : interpolation interval : max time
res_ts1 = resample(ts1, time, 'linear'); % Now resample it

A1=res_ts1.time; % new interpolated time
B1=res_ts1.data; % new interpolated parameter
tu = A1; % time
resamp = B1;       % interpolated dataset
variance = std(resamp)^2;
st_norm = (resamp - mean(resamp))/sqrt(variance) ;
% stdst = std([st_norm],1);
st_smooth = smooth(A1,st_norm,0.2,'loess'); %smooth the data
fn_resamp = st_norm - st_smooth; % Remove the trend (detrending)

%% Figure 1
figure (1)
set(gcf, 'Position',  [0, 0, 1000, 600])
subplot(311)
plot(age,parameter,'g',age,st1,'b'); xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(a) Parameter')
xticks(0:5:150)
hold off

subplot(312)
plot(tu,st_norm,tu,st_smooth); xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(b) Detrended parameter and trend')
xticks(0:5:150)
hold off

subplot(313)
plot(tu,fn_resamp); xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(c)Detrended parameter')
xticks(0:5:150)
hold off
% %% Requires R2020a or later
% exportgraphics(gcf,'vectorfig1.pdf','ContentType','vector')

%% Computation
%% Windowed
dt = tu(2)-tu(1) ;

% Calculate coi:
coi=coi_calc(fn_resamp,dt);

% Caulculate SST:
[Tff,w]=wsst(fn_resamp, 1./dt);

%% Inverse SST Individual Components:

freqrange1=[0.04, 0.0625]; % looking for 22 kyr period 1/22=0.0455
xrec1=iwsst(Tff,w,freqrange1); % blue line in the figure

freqrange2=[0.0625,0.125]; % looking for 11 kyr period 1/11=0.0909
xrec2=iwsst(Tff,w,freqrange2); % blue line in the figure

freqrange3=[0.1429,0.25]; % looking for 5.5 period 1/2.3=0.182
xrec3=iwsst(Tff,w,freqrange3); % green line in the figure

freqrange4=[0.3,0.56]; % looking for 2.3 kyr periods
xrec4=iwsst(Tff,w,freqrange4); % caena line in the figure

freqrange5=[0.6,1.25]; % looking for 1 kyr period 1/1=1
xrec5=iwsst(Tff,w,freqrange5); % magenda line in the figure

xrec=xrec1+xrec2+xrec3+xrec4+xrec5;

%% SST results for windowed frequencies
figure (2)
set(gcf, 'Position',  [0, 0, 1850, 600])

subplot(422) % manuscript Fig(c)
plot(tu,xrec1,'b'); hold on
plot(tu,xrec5,'m'); xlim([0 150])
mode1_v=max(xcorr(fn_resamp,xrec1)); 
mode5_v=max(xcorr(fn_resamp,xrec5)); 
%title(['Mode 1 reproduces ' num2str(floor(mode_v(1))) '% of Amplitude'])
axis([0 150 min(xrec5) max(xrec5)]); 
xticks(0:5:150); set(gca,'XTicklabels',[]); set(gca,'fontsize',12);
grid on; hold off;

subplot(424) % Fig(d)
plot(tu,xrec2,'r'); xlim([0 150])
mode2_v=max(xcorr(fn_resamp,xrec2)); 
%title(['Mode 2 reproduces ' num2str(floor(mode_v(2))) '% of Amplitude'])
xticks(0:5:150); set(gca,'XTicklabels',[]); set(gca,'fontsize',12);
grid on

subplot(426) % Fig(e)
plot(tu,xrec3,'g'); xlim([0 150])
mode3_v=max(xcorr(fn_resamp,xrec3)); 
%title(['Mode 3 reproduces ' num2str(floor(mode_v(3))) '% of Amplitude'])
% axis([0 150 -1.5 1.5]); 
xticks(0:5:150); set(gca,'XTicklabels',[]); set(gca,'fontsize',12);
grid on

subplot(428) % Fig(f)
plot(tu,xrec4,'c'); xlim([0 150])
mode4_v=max(xcorr(fn_resamp,xrec4)); 
%title(['Mode 4 reproduces ' num2str(floor(mode_v(4))) '% of Amplitude'])
axis([0 150 min(xrec4) max(xrec4)]); 
xticks(0:5:150); set(gca,'fontsize',12); xlabel('Time (Ka)')
grid on

total_var=0.28*floor(max(fn_resamp).*(mode1_v+mode2_v+mode3_v+mode4_v+mode5_v));

subplot(421) % Fig(a)
plot(tu,fn_resamp,tu,xrec); title(['Reproduces ' num2str(total_var) '% of Amplitude'])
xticks(0:5:150); 
axis([0 150 min(fn_resamp) max(fn_resamp)]); set(gca,'XTicklabels',[]); set(gca,'fontsize',12);
grid on

%% Noise Red Injection:
sigma = 0.25;%sqrt( var(clean)*10.^( -snrdb /10 ) );
% Sigma = 0.25 corresponds to 25% noise of 
% a signal to noise ratio of 4:1.
realization=100;

for i=1:realization;
    noise=rednoise(length(tu));
    stack=fn_resamp+(sigma.*noise);
    [Tf,w]=wsst(stack, 1./dt);    
    s_sgs(:,:,i)=Tf;
end
    
s_conf=sum(s_sgs,3)./realization;

%% Plot the Results
subplot(4,2,[3,5,7]) 
imagesc(tu, -log2(w), abs(Tff)); % manuscript Fig.4b
colormap (flipud(bone)); brighten(-0.9);
set(gca,'YDir','reverse'); ylabel('Period (Kyr)'); xlabel('Time (Ka)'); 
% set(gca,'XTicklabels',[]); 
yt=get(gca,'YTick'); yticks=2.^yt'; yticklabels({num2str(yticks)});
axis([0 150 -1 5]);
xticks(0:5:150); set(gca,'fontsize',12); xlabel('Time (Ka)');
% xticklabels({'0','','','','20','','','','40','','','','60','','','','80','','','','100','','','','120'})
% caxis([0 .05])
hold on
area(tu,log2(coi),8); alpha(0.5);
hold on
contour(tu,-log2(w),abs(s_conf),[-99,.002],'k');
hold on

%% Inverse SST search areas
yline(log2(1./freqrange1(1)),'b','LineWidth',2)
yline(log2(1./freqrange1(2)),'b','LineWidth',2)
yline(log2(1./freqrange2(1)),'r','LineWidth',2)
yline(log2(1./freqrange2(2)),'r','LineWidth',2)
yline(log2(1./freqrange3(1)),'g','LineWidth',2)
yline(log2(1./freqrange3(2)),'g','LineWidth',2)
yline(log2(1./freqrange4(1)),'c','LineWidth',2)
yline(log2(1./freqrange4(2)),'c','LineWidth',2)
yline(log2(1./freqrange5(1)),'m','LineWidth',2)
yline(log2(1./freqrange5(2)),'m','LineWidth',2)

%% Show periods as horizontal lines
yline(log2(1.0),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(1.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(2.75),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(5.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(11),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(22),'color',[0.75, 0.6, 0.1],'LineWidth',2)
hold off

%% Requires R2020a or later
% exportgraphics(gcf,'vectorfig2.pdf','ContentType','vector')

%%
%% Ridge Extraction not windowed

% Calculate coi
coi=coi_calc(fn_resamp,dt);

% Calculate SST:
[Tff,w] = wsst(fn_resamp, 1./dt);

% Ridges
[fridge,iridge] = wsstridge(Tff,6,w,'NumRidges',6);

% Inverse SST:
xrecc = iwsst(Tff,iridge); % find 

%% Plot the results for not windowed ridge extraction
figure (3)
set(gcf, 'Position',  [0, 0, 1850, 600])

subplot(4,2,[3,5,7])
imagesc(tu, -log2(w), abs(Tff)); % manuscript Fig.4h
colormap (flipud(bone)); brighten(-0.9);
set(gca,'YDir','reverse'); xlabel('Time (Ka)'); ylabel('Period (Myr)');
yt=get(gca,'YTick'); yticks=2.^yt'; yticklabels({num2str(yticks)});
axis([0 150 -1 5]);
xticks(0:5:150); set(gca,'fontsize',12); xlabel('Time (Ka)')
% caxis([0 .05])
hold on
area(tu,log2(coi),8); alpha(0.5);
hold on
plot(tu,-log2(fridge),'r--','linewidth',2); % marking ridges in Fig.4h

%% Show periods as horizontal lines
yline(log2(1.0),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(1.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(2.75),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(5.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(11),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(22),'color',[0.75, 0.6, 0.1],'LineWidth',2)
hold off

subplot(4,2,2) % manuscript Fig.4i
mode1=1;
plot(tu,xrecc(:,mode1)); xlim([0 150])
mode1_v=max(xcorr(fn_resamp,xrecc(:,mode1)));
%title(['Mode 1 reproduces ' num2str(floor(mode1_v*100)) '% of Variability'])
axis([0 150 min(xrecc(:,mode1)) max(xrecc(:,mode1))]); xticks(0:5:150);
set(gca,'fontsize',12); set(gca,'XTicklabels',[]);
grid on

subplot(4,2,4) % manuscript Fig.4j
mode3=3;
plot(tu,xrecc(:,mode3)); xlim([0 150])
mode2_v=max(xcorr(fn_resamp,xrecc(:,mode3)));
%title(['Mode 2 reproduces ' num2str(floor(mode2_v*100)) '% of Variability'])
axis([0 150 min(xrecc(:,mode3)) max(xrecc(:,mode3))]); xticks(0:5:150);
set(gca,'fontsize',12); set(gca,'XTicklabels',[]); 
grid on

subplot(4,2,6) % manuscript Fig.4k
mode5=5;
plot(tu,xrecc(:,mode5)); xlim([0 150])
mode5_v=max(xcorr(fn_resamp,xrecc(:,mode5)));
%title(['Mode 3 reproduces ' num2str(floor(mode3_v*100)) '% of Variability'])
axis([0 150 min(xrecc(:,mode5)) max(xrecc(:,mode5))]); xticks(0:5:150);
set(gca,'fontsize',12); set(gca,'XTicklabels',[]); 
grid on

subplot(4,2,8) % manuscript Fig.4k
mode2=2;
mode4=4;
mode6=6;
plot(tu,xrecc(:,mode2),'m',tu,xrecc(:,mode4),'b',tu,xrecc(:,mode6),'g'); xlim([0 150])
mode2_v=max(xcorr(fn_resamp,xrecc(:,mode2)));
mode4_v=max(xcorr(fn_resamp,xrecc(:,mode4)));
mode6_v=max(xcorr(fn_resamp,xrecc(:,mode6)));
%title(['Mode 3 reproduces ' num2str(floor(mode3_v*100)) '% of Variability'])
axis([0 150 min(xrecc(:,mode4)) max(xrecc(:,mode4))]); 
xticks(0:5:150); set(gca,'fontsize',12); xlabel('Time (Ka)')
grid on

total_var=0.15*floor(1*max(fn_resamp).*1*(mode1_v+mode2_v+mode3_v+mode4_v+mode5_v+mode6_v));

subplot(4,2,1) % addition to Fig.4g
plot(tu,fn_resamp,tu,sum(xrecc,2),'g'); 
axis([0 150 min(fn_resamp) max(fn_resamp)]); set(gca,'xticklabel',[])
title(['Reproduces ' num2str(total_var) '% of Amplitude'])
xticks(0:5:150); axis([0 150 -1.5 1.5]); set(gca,'fontsize',12);
grid on

%% Requires R2020a or later
% exportgraphics(gcf,'vectorfig3.pdf','ContentType','vector')

%%
%% Ridge Extraction windowed

% Calculate coi
coi = coi_calc(fn_resamp, dt);

% Calculate SST
[Tff, w] = wsst(fn_resamp, 1./dt);

% Define frequency ranges
freqranges = [
    0.04, 0.0625;     % Example for 23 kyr period
    0.0625, 0.125;  % Example for 11 kyr period
    0.1429, 0.25;     % Example for 5.5 kyr period
    0.3, 0.56;     % Example for 2.3 kyr period
    0.6, 1.25        % Example for 1 kyr period
];

% Initialize cell arrays to store filtered frequencies and ridges
xrecc = cell(size(freqranges, 1), 1);
strongest_ridges = cell(size(freqranges, 1), 1);

% Calculate filtered frequencies for each frequency range separately
for k = 1:size(freqranges, 1)
    % Find the indices corresponding to the frequency range
    freq_idx = find(w >= freqranges(k,1) & w <= freqranges(k,2));
    
    % Filter the time-frequency representation within the frequency range
    Tff_filtered = Tff(freq_idx, :);
    w_filtered = w(freq_idx);
    
    % Extract ridges within the filtered frequency range
    [fridge, iridge] = wsstridge(Tff_filtered, 5, w_filtered, 'NumRidges', 5);
    
    % Find the ridge with the maximum amplitude
    [~, max_ridge_idx] = max(max(abs(fridge), [], 1));
    strongest_ridge = fridge(:, max_ridge_idx);
    strongest_ridges{k} = strongest_ridge;
    
    % Reconstruct the signal using the ridge indices
    xrecc{k} = iwsst(Tff_filtered, iridge(:, max_ridge_idx));
end

%% Plot the results for windowed ridge extraction

figure(4)
set(gcf, 'Position', [0, 0, 1850, 600])

subplot(4,2,[3,5,7])
imagesc(tu, -log2(w), abs(Tff)); % manuscript Fig.4h
colormap(flipud(bone));
brighten(-0.9);
set(gca,'YDir','reverse');
xlabel('Time (Ka)');
ylabel('Period (Myr)');
yt = get(gca,'YTick');
yticks = 2 .^ yt';
yticklabels({num2str(yticks)});
axis([0 150 -1 5]);
xticks(0:5:150);
set(gca,'fontsize',12);
xlabel('Time (Ka)');
hold on
area(tu,log2(coi),8); % Highlight the Cone of Influence
alpha(0.5);

% Plot the strongest ridge for each frequency range
colors = {'b--', 'r--', 'g--', 'c--','m--'};
for k = 1:size(freqranges, 1)
    strongest_ridge = strongest_ridges{k};
    plot(tu, -log2(strongest_ridge), colors{k}, 'LineWidth', 2); % Plot strongest ridges
end
hold off

%% Show periods as horizontal lines
yline(log2(1.0),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(1.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(2.75),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(5.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(11),'color',[0.75, 0.6, 0.1],'LineWidth',2)
yline(log2(22),'color',[0.75, 0.6, 0.1],'LineWidth',2)

%% Plot straight lines for each frequency range
for k = 1:size(freqranges, 1)
    line([tu(1), tu(end)], -log2([freqranges(k,1), freqranges(k,1)]), 'Color', colors{k}(1), 'LineStyle', '-', 'LineWidth', 1.5);
    line([tu(1), tu(end)], -log2([freqranges(k,2), freqranges(k,2)]), 'Color', colors{k}(1), 'LineStyle', '-', 'LineWidth', 1.5);
end
hold off

%% Plot reconstructed signals for each frequency range
subplot(4,2,2)
plot(tu,xrecc{1},'b'); hold on
plot(tu,xrecc{5},'m'); xlim([0 150])
mode1_v=max(xcorr(fn_resamp,xrecc{1}));
mode15_v=max(xcorr(fn_resamp,xrecc{5}));
%title(['Mode 1 reproduces ' num2str(floor(mode1_v*100)) '% of Variability'])
axis([0 150 min(xrecc{5}) max(xrecc{5})]); xticks(0:5:150);
set(gca,'fontsize',12); set(gca,'XTicklabels',[]);
grid on; hold off

subplot(4,2,4) % 
plot(tu,xrecc{2},'r'); xlim([0 150])
mode2_v=max(xcorr(fn_resamp,xrecc{2}));
%title(['Mode 2 reproduces ' num2str(floor(mode2_v*100)) '% of Variability'])
axis([0 150 min(xrecc{2}) max(xrecc{2})]); xticks(0:5:150);
set(gca,'fontsize',12); set(gca,'XTicklabels',[]); 
grid on

subplot(4,2,6) %
plot(tu,xrecc{3},'g'); xlim([0 150])
mode3_v=max(xcorr(fn_resamp,xrecc{3}));
%title(['Mode 3 reproduces ' num2str(floor(mode3_v*100)) '% of Variability'])
axis([0 150 min(xrecc{3}) max(xrecc{3})]); xticks(0:5:150);
set(gca,'fontsize',12); set(gca,'XTicklabels',[]); 
grid on

subplot(4,2,8) % 
plot(tu,xrecc{4},'c'); xlim([0 150])
mode4_v=max(xcorr(fn_resamp,xrecc{4}));
%title(['Mode 3 reproduces ' num2str(floor(mode3_v*100)) '% of Variability'])
axis([0 150 min(xrecc{4}) max(xrecc{4})]); 
xticks(0:5:150); set(gca,'fontsize',12); xlabel('Time (Ka)')
grid on; hold off

%% Reconstructed signal power
total_var=0.3*floor(1*max(fn_resamp).*1*(mode1_v+mode2_v+mode3_v+mode4_v+mode5_v));

%% Plot the original and reconstructed signals
subplot(4,2,1) % addition to Fig.4g
plot(tu, fn_resamp, tu, sum(cat(2, xrecc{:}), 2), 'g'); 
axis([0 150 min(fn_resamp) max(fn_resamp)]); set(gca,'xticklabel',[])
title(['Reproduces ' num2str(total_var) '% of Amplitude'])
xticks(0:5:150); axis([0 150 -1.5 1.5]); set(gca,'fontsize',12);
grid on

%% Requires R2020a or later
% exportgraphics(gcf,'vectorfig4.pdf','ContentType','vector')

%% Plot SST results for windowed frequencies of ~11kyr, ~5.5kyr, ~2.75kyr
figure (5)
set(gcf, 'Position',  [0, 0, 1000, 600])

subplot(311)
plot(tu,fn_resamp,'color',[0.0, 0.0, 0.8]); hold on
plot(tu,2.5*xrec2,'r');
xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(a) 11 kyr cycle')
xticks(0:5:150)
hold off

subplot(312)
plot(tu,fn_resamp,'color',[0.0, 0.0, 0.8]); hold on
plot(tu,3.*xrec3,'g');
xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(b) 5.5 kyr cycle')
xticks(0:5:150)
hold off

subplot(313)
plot(tu,fn_resamp,'color',[0.0, 0.0, 0.8]); hold on
plot(tu,3.5*xrec4,'c'); hold on
xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(c) 2.3 kyr cycle')
xticks(0:5:150)
hold off

%% Requires R2020a or later
% exportgraphics(gcf,'vectorfig5.pdf','ContentType','vector')

%% Plot SST ridge extraction results for windowed frequencies of ~11kyr, ~5.5kyr, ~2.75kyr
figure (6)
set(gcf, 'Position',  [0, 0, 1000, 600])

subplot(311)
plot(tu,fn_resamp,'color',[0.0, 0.0, 0.8]); hold on
plot(tu,2.5*xrecc{2},'r');
xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(a) 11 kyr cycle')
xticks(0:5:150)
hold off

subplot(312)
plot(tu,fn_resamp,'color',[0.0, 0.0, 0.8]); hold on
plot(tu,3*xrecc{3},'g');
xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(b) 5.5 kyr cycle')
xticks(0:5:150)
hold off

subplot(313)
plot(tu,fn_resamp,'color',[0.0, 0.0, 0.8]); hold on
plot(tu,3.5*xrecc{4},'c');
xlim([0 150])
grid on
xlabel('Time (Ka)')
ylabel('Parameter (units)')
title('(c) 2.3 kyr cycle')
xticks(0:5:150)
hold off

%% Requires R2020a or later
% exportgraphics(gcf,'vectorfig6.pdf','ContentType','vector')

%% Hilbert Transform for Spectral Analysis
% Compute the Hilbert Transform to analyze signal phase and amplitude
%% Initialize Hilbert Spectrum
H = [];
% Apply Hilbert transform to each "xrec" component
xrec_all = {xrec1, xrec2, xrec3, xrec4, xrec5}; % Store components in a cell array for looping

for i = 1:length(xrec_all)
    current_component = xrec_all{i}; % Extract current component
    
    % Hilbert transform of the i-th component
    ht = hilbert(current_component);
    
    % Instantaneous amplitude and frequency
    amplitude = abs(ht);
    phase = unwrap(angle(ht));
    frequency = diff(phase) / (2 * pi * mean(diff(tu)));  % Calculate instantaneous frequency
    
    % Resample frequency to match the length of the time vector
    frequency = [frequency; frequency(end)];
    
    % Append to Hilbert Spectrum
    H = [H; amplitude .* frequency];
end

%% Plot the Hilbert Spectrum
figure (7);
imagesc(tu, 1:size(H, 1), H);
xlabel('Time (Ka)');
ylabel('Component Index');
title('Hilbert Spectrum of Synchrosqueezing Components');
colorbar;
set(gca, 'YDir', 'normal');

%% Corrected reconstruction of the signal from the xrec components
reconstructed_signal = zeros(size(tu)); % Initialize reconstructed signal

% Loop through the xrec components and sum them
for i = 1:length(xrec_all)
    reconstructed_signal = reconstructed_signal + xrec_all{i}; % Sum each component
end

%% Plot the original signal and the sum of the reconstructed components for reference
figure (8);
set(gcf, 'Position',  [0, 0, 1000, 600])
subplot(4,1,1);
plot(tu, fn_resamp); hold on;
plot(tu, reconstructed_signal);
legend('Original Signal', 'Reconstructed Signal from xrec Components');
xlabel('Time (Ka)');
ylabel('Amplitude');
title('Original and Reconstructed Signal');
xticks(0:5:150)
grid on;

%% Plot the individual components
figure (9);
set(gcf, 'Position',  [0, 0, 1000, 600])
for i = 1:length(xrec_all)
    subplot(length(xrec_all), 1, i);
    plot(tu, xrec_all{i});
    ylabel(['Component ', num2str(i)]);
    xticks(0:5:150)
    grid on;
end
xlabel('Time (Ka)');
sgtitle('Reconstructed Components from Synchrosqueezing Transform');

%% Parameters for significance testing
alpha = 0.05;  % Significance level
n = length(tu);  % Length of the time vector (tu)

% Loop through each xrec component for power spectrum analysis
for i = 1:length(xrec_all)
    current_component = xrec_all{i}; % Extract current xrec component
    
    % Calculate periodogram
    [Pxx, F] = periodogram(current_component, [], [], 1/mean(diff(tu)));  % Periodogram of current component
    
    % White noise theoretical spectrum (expected value)
    dof = 2;  % Degrees of freedom for periodogram (for each component)
    chi2crit = chi2inv(1 - alpha, dof);  % Chi-square critical value for the significance level
    sig95 = chi2crit * mean(Pxx) / dof;  % 95% significance level
    
    %% Power/Frequency Hilbert Transform Plot (in dB scale)
    figure(11 + (i-1));  % Create a new figure for each component
    set(gcf, 'Position',  [0, 0, 1000, 600])
    % Plot the power spectrum in dB scale
    subplot(211)
    plot(F, 10*log10(Pxx), 'b', 'LineWidth', 1.5);  % Convert Power to dB scale
    hold on;
    % Plot the 95% significance level line (Power/Frequency in dB)
    plot(F, 10*log10(repmat(sig95, size(F))), 'r--', 'LineWidth', 1.5);
    axis([0 2 min(10*log10(Pxx)) max(10*log10(Pxx))]);  % Adjust the axis limits as needed
    xticks(0:0.1:2);
    xline((1/22),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/11),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/5.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/2.75),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/1),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xlabel('Frequency (1/kyr)');
    ylabel('Power/Frequency (dB/Hz)');
    title(['Power Spectrum (Power/Frequency) of Component ', num2str(i)]);
    grid on;
    
    % Add a legend
    legend('Power Spectrum', 'Significance Level (sig95)');
    
    %% Hilbert Transform FFT-like Power-Only Plot
    subplot(212)
    % figure(21 + (i-1));  % Create a separate figure for the Power plot (not in dB)
    % Plot the power spectrum (Power only)
    plot(F, Pxx, 'b', 'LineWidth', 1.5);
    hold on;
    % Plot the significance level line (Power only)
    plot(F, repmat(sig95, size(F)), 'r--', 'LineWidth', 1.5);
    axis([0 2 min(Pxx) max(Pxx)]);  % Adjust the axis limits as needed
    xticks(0:0.1:2);
    xline((1/22),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/11),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/5.5),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/2.75),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xline((1/1),'color',[0.75, 0.6, 0.1],'LineWidth',2)
    xlabel('Frequency (1/kyr)');
    ylabel('Power');
    title(['Power Spectrum (Power) of Component ', num2str(i)]);
    grid on;
    % Add a legend
    legend('Power Spectrum', 'Significance Level (sig95)');
end


%% The end of the code
