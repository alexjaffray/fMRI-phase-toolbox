function logElapsed(cfg, label, timerValue)
%LOGELAPSED Print elapsed time for a pipeline step.

if isfield(cfg, 'verbose') && cfg.verbose
    timestamp = string(datetime("now", "Format", "HH:mm:ss"));
    fprintf('[%s] %s finished in %.1f s.\n', ...
        timestamp, label, toc(timerValue));
end

end
