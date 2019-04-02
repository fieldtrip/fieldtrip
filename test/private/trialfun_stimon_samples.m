function trl = trialfun_stimon_samples(cfg)
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
correctresponse  = 10041;
begintrial       = 10044;
endtrial         = 10045;
stimon           = 10030;
distractorChange = 12000;
targetChange     = 12001;
attCnds          = 20001:20004; % att in/out by target change first/second
E          = struct2cell(event);
samples    = cell2mat(E(1,:)); % now in vector form
value      = cell2mat(E(2,:));
begmark    = find(value==begintrial); % loop through the trial beginnings
endmark    = find(value==endtrial); % loop through the trial beginnings
trl        = []; % initialize the cfg.trl
for k=1:length(begmark)
    vals = value(begmark(k):endmark(k));
    if any(ismember(vals,attCnds)) && ~isempty(find(vals==correctresponse))
        % create the trl matrix in sample units
        samp = samples(begmark(k):endmark(k)); % in timestamp units	
        beginSamp      = samp(find(vals==stimon));        
        sampDistractor = samp(find(vals==distractorChange));
        sampTarget     = samp(find(vals==targetChange));       
        endSamp        = min([sampTarget(:);sampDistractor(:)]); % limit until first change        
        offset         = -round(hdr.Fs*2.75);        
        trl            = [trl; [beginSamp+offset endSamp offset]];
    end 
end
