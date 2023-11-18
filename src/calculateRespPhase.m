function respPhase = calculateRespPhase(pulset, time, doHilbert)

% adapted from TAPAS

verbose.level = 2;

resp_max = inf;

overshoot = find(abs(pulset) > resp_max);
pulset(overshoot) = resp_max;

rsampint = time(2) - time(1);
% Compute normalized signal and the sign of the derivative

% pulset = tapas_physio_filter_respiratory(pulset, rsampint, [0.01, 0.869], false, false, verbose);

maxr = max(pulset);
minr = min(pulset);

npulse = (pulset-minr)/(maxr-minr);

% Calculate derivative of normalised pulse wrt time
% over 1 sec of data as described in Glover et al.
ksize = round(0.5 * (1/rsampint));
kernel = ksize:-1:-ksize; kernel = kernel ./ sum(kernel.^2);
dpulse = tapas_physio_conv(pulset, kernel, 'symmetric');

% This uses a quadratic Savitzky-Golay filter, for which the coefficients
% have a simple linear form. See e.g.
% https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
% `scipy.signal.savgol_coeffs(5, polyorder=2, deriv=1, use='conv')`

% Tolerance to the derivative
% dpulse(abs(dpulse) < 1e-6) = 0;
dpulse = sign(dpulse);

% number of histogram bins determined by dynamic range of detected values
% and length of input time course
nbins = min(length(unique(pulset)), floor(length(pulset)/100));

[h, rout] = hist(npulse(dpulse ~=0 & npulse < resp_max), nbins);

binnum = floor(npulse*(nbins-1)) + 1;
binnum(overshoot) = nbins;

cumsumh = cumsum(h');
sumh = cumsumh(end);

rphase = pi*(cumsumh(binnum)/sumh).*dpulse+pi;

respPhase = rphase;

normalized_trace = npulse - 0.5;

if doHilbert
    
    t = time; % Time vector (seconds)
    fs = 1/(t(2) - t(1));
    respiratory_signal = normalized_trace;
    
%     % Plot the respiratory signal
%     figure;
%     subplot(2,1,1);
%     plot(t, respiratory_signal);
%     title('Respiratory Signal');
%     xlabel('Time (s)');
%     ylabel('Amplitude');
    
    % Apply hilbert respiratory phase estimation method
    hilbert_transform = hilbert(respiratory_signal);
    instantaneous_phase = angle(hilbert_transform);
    
    % Convert phase to the range [0, 2*pi]
    instantaneous_phase = mod(instantaneous_phase, 2*pi);
%     
%     % Plot the instantaneous phase
%     subplot(2,1,2);
%     plot(t, instantaneous_phase);
%     title('Instantaneous Phase');
%     xlabel('Time (s)');
%     ylabel('Phase (rad)');
    
    % Optionally, extract the respiratory frequency (we don't do this yet)
    instantaneous_frequency = diff(unwrap(instantaneous_phase)) * fs / (2*pi);
%     figure;
%     plot(t(1:end-1), instantaneous_frequency);
%     title('Instantaneous Respiratory Frequency');
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
    

respPhase = instantaneous_phase;
    
end


end



