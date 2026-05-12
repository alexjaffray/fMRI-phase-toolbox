function phase = estimatePhase(signal, time, opts)
%ESTIMATEPHASE Estimate respiratory phase from a respiratory timecourse.
%
% The default Hilbert method returns phase in radians on [0, 2*pi).

arguments
    signal (:,1) double
    time (:,1) double
    opts.method (1,1) string = "hilbert"
end

if numel(signal) ~= numel(time)
    error('fmriPhase:respiration:PhaseLengthMismatch', ...
        'signal and time must have the same number of samples.');
end

signal = signal(:);

if numel(signal) < 3
    phase = nan(size(signal));
    return
end

if all(~isfinite(signal)) || range(signal, 'omitnan') == 0
    phase = zeros(size(signal));
    return
end

signal = fillmissing(signal, 'linear', 'EndValues', 'nearest');
signal = signal - mean(signal, 'omitnan');

switch lower(opts.method)
    case "hilbert"
        analyticSignal = hilbert(signal);
        phase = mod(angle(analyticSignal), 2 * pi);
    otherwise
        error('fmriPhase:respiration:UnsupportedPhaseMethod', ...
            'Unsupported respiratory phase method "%s".', opts.method);
end

end
