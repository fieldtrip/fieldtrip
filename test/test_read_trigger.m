function test_read_trigger

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_event read_trigger

hdr = [];
hdr.Fs = 1000;
hdr.nChans = 1;
hdr.nSamples = 100;
hdr.nSamplesPre = 0;
hdr.label = {'1'};

dat = zeros(hdr.nChans, hdr.nSamples);

filename = tempname;
matfile = [filename '.mat'];
binfile = [filename '.bin'];

%%

dat(:) = 0;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1);
assert(numel(event)==0);

%%

dat(:)  = 0;
dat(51) = 1;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'up');

assert(numel(event)==1);
assert(event.sample==51);

%%

dat(:) = 1;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'up', 'trigpadding', true);
assert(numel(event)==0); % it will not be detected in this case
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'up', 'trigpadding', false);
assert(numel(event)==1); % it will be detected in this case
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'down', 'trigpadding', true);
assert(numel(event)==0); % it will not be detected in this case
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'down', 'trigpadding', false);
assert(numel(event)==1); % it will be detected in this case

%%

dat(:)  = 0;
dat(51) = 1;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'down');

assert(numel(event)==1);
assert(event.sample==52);

%%

dat(:)  = 0;
dat(51) = 1;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'both');

assert(numel(event)==2);
assert(event(1).sample==51);
assert(event(2).sample==52);

%%

dat(:)  = 0;
dat(51) = 1;
dat(52) = 2;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'updiff');

% there are two going up
assert(numel(event)==2);
assert(event(1).sample==51);
assert(event(2).sample==52);

event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'downdiff');

% and one going down
assert(numel(event)==1);
assert(event(1).sample==53);

%%

dat(:)  = 0;
dat(51) = -1;
dat(52) = -2;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'downdiff');

% there are two going down
assert(numel(event)==2);
assert(event(1).sample==51);
assert(event(2).sample==52);

event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'updiff');

% and one going up
assert(numel(event)==1);
assert(event(1).sample==53);

%%

dat(:)  = 0;
dat(51) = 1;
dat(52) = 2;
dat(53) = 3;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'biton');

assert(numel(event)==3);
assert(event(1).sample==51 && event(1).value==1); % bit 1 goes on
assert(event(2).sample==52 && event(2).value==2); % bit 1 goes off, 2 goes on
assert(event(3).sample==53 && event(3).value==1); % bit 1 goes on

%%

dat(:)  = 0;
dat(51) = 1;
dat(52) = 2;
dat(53) = 3;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'bitoff');

assert(event(1).sample==52 && event(1).value==1); % bit 1 goes off
assert(event(2).sample==54 && event(2).value==1); % bit 1 goes off
assert(event(3).sample==54 && event(3).value==2); % bit 2 goes off as well

%%

dat(:) = 0;
dat(1:2:end) = 1;
dat(2:2:end) = 0;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'up', 'denoise', false);
assert(numel(event)==49); % the first one is not detected due to trigpadding
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'down', 'denoise', false);
assert(numel(event)==50);

%%

dat(:) = 0;
dat(1:2:end) = 0;
dat(2:2:end) = 1;
ft_write_data(matfile, dat, 'header', hdr, 'dataformat', 'fcdc_matbin');
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'up', 'denoise', false);
assert(numel(event)==50);
event = ft_read_event(matfile, 'chanindx', 1, 'detectflank', 'down', 'denoise', false);
assert(numel(event)==49); % the last one is not detected due to trigpadding

%%

delete(matfile)
delete(binfile)
