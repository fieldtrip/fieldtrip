function test_bug3368

% WALLTIME 00:10:00
% MEM 4gb
% DEPENDENCY

% load a fieldtrip dataset

ieeg_name      = dccnpath(fullfile('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog','SubjectNY394','NY394_VisualLoc_R1.edf'));
cfg            = [];
cfg.dataset    = ieeg_name;
cfg.continuous = 'yes';
cfg.channel    = 'all';
ft_data        = ft_preprocessing(cfg);

% name to write the file
ieeg_name_gdf = [tempname, '.gdf'];

% fetch the header
hdr_data = ft_fetch_header(ft_data);

% write the data
ft_write_data(ieeg_name_gdf,ft_data.trial{1},'header',hdr_data,'dataformat','gdf')

% load the data back in
test_ft_data      = ft_read_data(ieeg_name_gdf);
test_ft_header    = ft_read_header(ieeg_name_gdf);

% now check whether the loaded data are the same as the written data,
% the following should be zero:
[ok,msg] = isalmostequal(ft_data.trial{1}, test_ft_data);

indx = find(isfinite(test_ft_data));
assert(isequal(test_ft_data(indx), ft_data.trial{1}(indx))); % ensure that at least all numerical data is equal
