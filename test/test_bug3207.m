function test_bug3207

% WALLTIME 00:20:00
% MEM 4gb

% TEST ft_read_event read_edf

filename = 'test3_2048Hz.EDF';


hdr1 = ft_read_header(filename, 'chanindx', 1);
dat1 = ft_read_data(filename, 'header', hdr1);
tim1 = (1:hdr1.nSamples)/hdr1.Fs;


hdr2 = ft_read_header(filename, 'chanindx', 2);
dat2 = ft_read_data(filename, 'header', hdr2);
tim2 = (1:hdr2.nSamples)/hdr2.Fs;



evt1 = ft_read_event(filename, 'detectflank', []); % expressed at 45 Hz
evt2 = evt1;
for i=1:numel(evt2)
  evt2(i).sample = evt2(i).timestamp*hdr2.Fs + 1;  % expressed at 2048 Hz
end

figure
hold on
plot(tim2, dat2 ./ max(dat2), 'g');
plot([evt1.timestamp], ones(1,numel(evt1)), 'ro');
plot([evt2.timestamp], ones(1,numel(evt2)), 'm+');

