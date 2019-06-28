function cfg = save_large_cfg_fields(cfg, highest_ft, reproducescript_dir, nowtime)

% SAVE_LARGE_CFG_FIELDS is a helper function for ft_postamble_savevar and ft_postamble_savefig, and
% is used for the cfg.reproducescript functionality.

fn = ignorefields('recursesize');
for i=1:numel(fn)
  if isfield(cfg.callinfo.usercfg, fn{i}) && isstruct((cfg.callinfo.usercfg.(fn{i})))
    outputfile = make_or_fetch_inputfile(reproducescript_dir,...
      sprintf('%s_%s_input_%s.mat', nowtime, highest_ft, ...
      fn{i}), fn{i}, cfg.callinfo.usercfg.(fn{i}));
    cfg.callinfo.usercfg.(fn{i})  = outputfile;
  end
end

end