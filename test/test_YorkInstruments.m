function test_YorkInstruments

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_filetype ft_read_header ft_read_data ft_read_event

%%
% this is where the data is installed on the DCCN compute cluster
% and where the nightly regression testing is done

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/YI.meghdf5');
filename = '/mnt/megdata-yi/E3002/20210429_C1023/processed/20210429130007.meghdf5'
hdr = ft_read_header(filename);

if length(ft_channelselection('MEG', hdr.label))~=246
  error('did not select all MEG channels');
end

if ~ft_senstype(hdr, 'bti248')
  error('hdr ~= bti248');
end

if ~strcmp(ft_filetype(filename), 'York_Instruments_hdf5')
  error('Filetype not detected correctly');
end	

