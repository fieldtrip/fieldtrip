function inspect_issue1674

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_prepare_layout ft_convert_coordsys

%%

elec = [];
elec.label = {
  'Fpz'
  'T7'
  'Cz'
  'T8'
  'Oz'
  };
elec.elecpos = [
   80   0   0 % a
    0  80   0 % l
    0   0  80 % s
    0 -80   0
  -80   0   0
  ];
elec.unit = 'mm';
elec.coordsys = 'eeglab';

%%

% this was the key issue that was reported, Fpz should be at the front of the head

cfg = [];
cfg.elec = elec;
layout = ft_prepare_layout(cfg);
ft_plot_layout(layout);


%%

% this was a side issue that I detected while testing, since 'eeglab' is also 'als'
% the conversion from 'als' to 'ras' (and others) was not working as expected 

elec = [];
elec.label = {
  'a'
  'l'
  's'
  };
elec.elecpos = [
  1 0 0 % a
  0 1 0 % l
  0 0 1 % s
  ];
elec.unit = 'mm';
elec.coordsys = 'als';

target = {'als' 'ali' 'ars' 'ari' 'pls' 'pli' 'prs' 'pri' 'las' 'lai' 'ras' 'rai' 'lps' 'lpi' 'rps' 'rpi' 'asl' 'ail' 'asr' 'air' 'psl' 'pil' 'psr' 'pir' 'sal' 'ial' 'sar' 'iar' 'spl' 'ipl' 'spr' 'ipr' 'sla' 'ila' 'sra' 'ira' 'slp' 'ilp' 'srp' 'irp' 'lsa' 'lia' 'rsa' 'ria' 'lsp' 'lip' 'rsp' 'rip'};

elec = {elec};

for i=2:length(target)
  disp(target{i});
  
  elec{i} = ft_convert_coordsys(elec{i-1}, target{i});
  
  % the 1st electrode is at the nose
  letter1indx = find(elec{i}.elecpos(1,:)~=0); % where did it go?
  if elec{i}.elecpos(1,letter1indx)==1         % is the axis flipped?
    letter1val = 'a';
  else
    letter1val = 'p';
  end
  % the 2nd electrode is at the left
  letter2indx = find(elec{i}.elecpos(2,:)~=0);
  if elec{i}.elecpos(2,letter2indx)==1
    letter2val = 'l';
  else
    letter2val = 'r';
  end
  % the 3rd electrode is at the top
  letter3indx = find(elec{i}.elecpos(3,:)~=0);
  if elec{i}.elecpos(3,letter3indx)==1
    letter3val = 's';
  else
    letter3val = 'i';
  end
  
  result = 'xxx';
  result(letter1indx) = letter1val;
  result(letter2indx) = letter2val;
  result(letter3indx) = letter3val;
  
  assert(isequal(result, elec{i}.coordsys));
end



