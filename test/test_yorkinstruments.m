function test_yorkinstruments

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_filetype ft_read_header ft_read_data ft_read_event

%%
% this is where the data is installed on the DCCN compute cluster
% and where the nightly regression testing is done

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/yorkinstruments/yi.meghdf5');
hdr = ft_read_header(filename);

if length(ft_channelselection('MEG', hdr.label))~=246
  error('did not select all MEG channels');
end

if ~ft_senstype(hdr, 'bti248')
  error('hdr ~= bti248');
end

if ~strcmp(ft_filetype(filename), 'yorkinstruments_hdf5')
  error('Filetype not detected correctly');
end	

