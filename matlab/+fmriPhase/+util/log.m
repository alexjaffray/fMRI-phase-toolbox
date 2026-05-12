function log(cfg, message)
%LOG Print a timestamped pipeline message when cfg.verbose is true.

if isfield(cfg, 'verbose') && cfg.verbose
    timestamp = string(datetime("now", "Format", "HH:mm:ss"));
    fprintf('[%s] %s\n', timestamp, message);
end

end
