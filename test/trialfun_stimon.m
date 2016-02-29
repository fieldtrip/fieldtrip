function trl = trialfun_stimon(cfg)
 
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
correctresponse  = 10041;
begintrial       = 10044;
endtrial         = 10045;
stimon           = 10030;
distractorChange = 12000;
targetChange     = 12001;
attCnds          = 20001:20004;       % att in/out by target change first/second
E          = struct2cell(event);
samples    = cell2mat(E(1,:));        % now in vector form
value      = cell2mat(E(2,:));
timestamps = cell2mat(E(3,:));        % now in vector form
begmark    = find(value==begintrial); % loop through the trial beginnings
endmark    = find(value==endtrial);   % loop through the trial beginnings
trl = [];
for k=1:length(begmark)
  vals = value(begmark(k):endmark(k));
  if any(ismember(vals,attCnds)) && ~isempty(find(vals==correctresponse))
    ts = timestamps(begmark(k):endmark(k)); % in timestamp units
    beginTs      = ts(find(vals==stimon));
    tsDistractor = ts(find(vals==distractorChange));
    tsTarget     = ts(find(vals==targetChange));
    endTs        = min([tsTarget(:);tsDistractor(:)]);    % limit until first change
    offset       = - hdr.Fs*hdr.TimeStampPerSample*2.75;  % 40000 timestamps per second x 2.75 sec
    trl          = [trl; [beginTs+offset endTs offset]];
  end
end
