
%% Pipeline for getting breathing belt regressor for fMRI processing
%
% Use: Run the script, a dialog box will open and you can choose the
% physlog file you want to process. It will process and plot, and save the regressor matrix.  then you
% can load the breathing reg matrix into conn
% 15 May 2024
% Editors (alphabetical order):
% Alexander Jaffray: ajaffray@physics.ubc.ca

clear all;
close all;

%%
addpath('/srv/data/ajaffray/fMRI-phase-toolbox/');

%% Ask User for the phys log file!
[logFile,logDir] = uigetfile("*.log","Select the Relevant PhysLog File");

%%
physLogTable = readPhysLog(fullfile(logDir,logFile));

%%
load pipeline/scan_params.mat

%% Prepare physLog
pLogFileName = string(fullfile(logDir,logFile));
logfiles.cardiac = pLogFileName;
logfiles.respiration = pLogFileName;
logfiles.sampling_interval = 1/500;
logfiles.relative_start_acquisition = 0;
phaseSign = -1; % NEED
doTRComp = 0; % NEED

if phaseSign == -1
    doTRComp = 1;
end

[c, r, t2, cpulse, acq_codes] = tapas_physio_read_physlogfiles_philips(logfiles, 'PPU');

c = c(1:2:end);
r = r(1:2:end);
t2 = t2(1:2:end)/2;
scanStart = find(physLogTable.mark>3);
t2 = t2(scanStart+1:end);
[physLogResp,fh] = tapas_physio_filter_respiratory(r(scanStart+1:end),t2(2)-t2(1),[],true,true); 

samplingRate = length(physLogResp) / scanTime;
physLogTime = linspace(0,scanTime,length(physLogResp)); % NEED

%% Hilbert transform respiratory phase or retroicor resp phase computed separately

ts_physlog_derived = timeseries(physLogResp, physLogTime);
dt = 1.15;

physlog_derived_times = physLogTime(1):dt:physLogTime(end-1);

ts_physlog_derived_rs = resample(ts_physlog_derived, physlog_derived_times);


respPhasePhysLog_hf = calculateRespPhase(ts_physlog_derived.Data, ts_physlog_derived.Time, true);

breath_reg.reg = ts_physlog_derived_rs.Data;
save('breath_reg.mat','-struct','breath_reg')

%%
figure();
plot(ts_physlog_derived.Time,respPhasePhysLog_hf - pi);
ylabel('Respiratory Phase [rad]');
ylim([-4.5, 4.5]);
ax = gca;
set(ax,'FontWeight','bold','Box','on','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','YGrid','off','LineWidth',1.5, 'FontSize',18);
set(ax.XAxis,'TickDirection','in');
set(ax.YAxis,'TickDirection','out');
xlim([0 310]);
xlabel('Scan Time [s]');
hold on;
% Load the physlog file and process as before (naive approach)
resampledPhysLogTrace = physLogResp; %resample(physLogResp,interpolationFactor,round(samplingRate*TR)); % resample to same rate as the fmri-derived resp data
filteredPhysLogTrace = lowpass(resampledPhysLogTrace,20,500); % low pass
% filter to remove weird data jumps
respPhasePhysLog = calculateRespPhase(ts_physlog_derived_rs.Data, ts_physlog_derived_rs.Time, true);
plot(ts_physlog_derived_rs.Time,respPhasePhysLog - pi,'LineStyle', '-.', 'LineWidth', 2, 'Color', 'r');

%%

breath_reg.reg = respPhase(1:260);
save('breath_reg.mat','-struct','breath_reg')