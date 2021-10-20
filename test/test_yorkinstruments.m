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
  error('filetype not detected correctly');
end

%%

% there is no dewar to head transform available
sens = ft_read_sens(filename, 'senstype', 'meg', 'coordsys', 'dewar');
figure
ft_plot_sens(sens);

%%

headshape = ft_read_headshape(filename);
figure
ft_plot_headshape(headshape);

