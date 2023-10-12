function test_fieldtrip2ctf

% MEM 1gb
% WALLTIME 00:10:00
% DATA public
% DEPENDENCY fieldtrip2ctf

% Use the 'Subject01.ds' dataset for basic testing of the writing and reading

inputfile = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds');

%%

cfg         = [];
cfg.dataset = inputfile;
raw_input = ft_preprocessing(cfg);

% ds = readCTFds(inputfile);
% ds = raw_input.hdr.orig;

%%

outputfile = [tempname '.ds'];

%%

fieldtrip2ctf(outputfile, raw_input, 'ds', [])

%%

cfg         = [];
cfg.dataset = outputfile;
raw_output = ft_preprocessing(cfg);

%%

trialsel = 1;
chansel = 1;

figure;
subplot(2,1,1); plot(raw_input.time{trialsel}, raw_input.trial{trialsel}(chansel,:));
subplot(2,1,2); plot(raw_output.time{trialsel}, raw_output.trial{trialsel}(chansel,:));

%%

for trialsel=1:266
  chansel = 1:187;
  isalmostequal(raw_input.trial{trialsel}(chansel,:), raw_output.trial{trialsel}(chansel,:));
end

%%

% these are not the same

evt_input = ft_read_event(inputfile);
evt_output = ft_read_event(outputfile);

%%

% clean up

rmdir(outputfile, 's');

