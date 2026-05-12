%% Comparison of SVD components between subjects
% Michelle Medina
% August 8, 2022

%% Load saved data

close all;
clear all;

sub_M03 = load('sub_M03_fl_var.mat');
sub_M08 = load('sub_M08_fl_var.mat');
sub_M29 = load('sub_M29_fl_var.mat');

%% Plot histograms by plotting the time-course of the 2nd SVD coefficient and identify minima

% identify the minima (for gating of breath for example)
foundMins_sub_M03 = islocalmin(sub_M03.comp2,'MinProminence',0.0015);
foundMins_sub_M08 = islocalmin(sub_M08.comp2,'MinProminence',0.0015);
foundMins_sub_M29 = islocalmin(sub_M29.comp2,'MinProminence',0.0015);

% Plot
figure();
plot(sub_M03.timeVector,sub_M03.comp2);
hold on;
scatter(sub_M03.timeVector(foundMins_sub_M03),sub_M03.comp2(foundMins_sub_M03),'o');
title('2nd SVD Coefficient');
xlabel('Time (s)');
ylabel('Fluctuation in X map (ppm)');
legend('Delta X','Minima');

%% Make a histogram of the time intervals

% Gets the differences between adjacent minima times
differences_M03 = diff(sub_M03.timeVector(foundMins_sub_M03));
differences_M08 = diff(sub_M08.timeVector(foundMins_sub_M08));
differences_M29 = diff(sub_M29.timeVector(foundMins_sub_M29));

% Plot histogram for the different subjects
figure();
histogram(differences_M03, 10, 'BinLimits', [0,30]); % add BinLimits
hold on;
figure();
histogram(differences_M08, 10, 'BinLimits', [0,30]); % add BinLimits
hold on;
figure();
histogram(differences_M29, 10, 'BinLimits', [0,30]); % add BinLimits

%% Find the frequency of the 2nd SVD component 

% Constants
dwellTime = 1.15; % found in nifti header, 4th "PixelDimension"
Fs = 1/dwellTime; % sampling frequency
intFact = 1; % used interpolation factor
L = 260*intFact; % Length of signal
nsteps = 251;

% Compute Fourier Transforms
d_FT_M03 = fft(sub_M03.comp2,nsteps);
d_FT_M08 = fft(sub_M08.comp2,nsteps);
d_FT_M29 = fft(sub_M29.comp2,nsteps);

% Power spectrum
d_power_M03 = d_FT_M03.*conj(d_FT_M03)/nsteps;
d_power_M08 = d_FT_M08.*conj(d_FT_M08)/nsteps;
d_power_M29 = d_FT_M29.*conj(d_FT_M29)/nsteps;

% Frequency domain
f = Fs*(0:(nsteps-1))/nsteps;

% Plot
plot(f, d_power_M03);
hold on;
plot(f, d_power_M08);
hold on;
plot(f, d_power_M29);

% Use Periodogram to find max peak
% Male subject 3
[Pxx_M03, omega_M03] = periodogram(sub_M03.comp2,rectwin(length(sub_M03.comp2)),nsteps,Fs);
[psdPeaks_M03, peakLocations_M03] = findpeaks(Pxx_M03,omega_M03);

% Male subject 8
[Pxx_M08, omega_M08] = periodogram(sub_M08.comp2,rectwin(length(sub_M08.comp2)),nsteps,Fs);
[psdPeaks_M08, peakLocations_M08] = findpeaks(Pxx_M08,omega_M08);

% Male subject 29
[Pxx_M29, omega_M29] = periodogram(sub_M29.comp2,rectwin(length(sub_M29.comp2)),nsteps,Fs);
[psdPeaks_M29, peakLocations_M29] = findpeaks(Pxx_M29,omega_M29);

% Print max peak value
[m3,i3] = max(psdPeaks_M03)
[m8,i8] = max(psdPeaks_M08)
[m29,i29] = max(psdPeaks_M29)

% Location of max peak value in frequency
peakLocations_M03(i3)
peakLocations_M08(i8)
peakLocations_M29(i29)
