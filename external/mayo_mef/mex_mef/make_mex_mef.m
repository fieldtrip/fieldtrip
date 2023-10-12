% Compile mex files required to process MEF files

% Copyright 2019-2020 Richard J. Cui. Created: Wed 05/29/2019  9:49:29.694 PM
% $Revision: 0.8 $  $Date: Fri 06/05/2020 10:27:44.467 AM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905, USA
%
% Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)

% =========================================================================
% find the directory of make_mex_mef file
% =========================================================================
% go to mex_mef
mex_mef = fileparts(mfilename('fullpath')); % directory of make_mex_mef.m assumed in mex_mef
mayo_mef = fileparts(mex_mef); % mayo_mef directory 

% =========================================================================
% processing mex for MEF 2.0 file
% =========================================================================

% =========================================================================
% processing mex for MEF 2.1 file
% =========================================================================
cd([mex_mef,filesep,'mef_2p1']) % assume 'mef_2p1' is the subdirectory

fprintf('===== Compiling c-mex for MEF 2.1 data =====\n')
fprintf('Building read_mef_header_2p1.mex*\n')
mex -output read_mef_header_2p1 ...
    read_mef_header_mex_2p1.c ...
    mef_lib_2p1.c
movefile('read_mef_header_2p1.mex*',mayo_mef)

fprintf('\n')
fprintf('Building decompress_mef_2p1.mex*\n')
mex -output decompress_mef_2p1 ...
    decompress_mef_mex_2p1.c ...
    mef_lib_2p1.c
movefile('decompress_mef_2p1.mex*',mayo_mef)

cd(mayo_mef)

% =========================================================================
% processing mex for MEF 3.0 file
% =========================================================================
fprintf('\n')
fprintf('===== Compiling c-mex for MEF 3.0 data =====\n')

% cd([mex_mef, filesep, 'mef_3p0'])
% fprintf('Building read_mef_header_mex_3p0.mex*\n')
% mex -output read_mef_header_mex_3p0 ...
%     read_mef_header_mex_3p0.c ...
%     matmef/matmef_mapping.c ...
%     matmef/mex_datahelper.c
% movefile('read_mef_header_mex_3p0.mex*',mayo_mef)

cd([mex_mef, filesep, 'mef_3p0', filesep, 'matmef', filesep]) % assume 'mef_3p0' is the subdirectory
fprintf('\n')
fprintf('Building read_mef_session_metadata.mex*\n')
mex -output read_mef_session_metadata ...
    read_mef_session_metadata.c ...
    matmef_mapping.c ...
    mex_datahelper.c
movefile('read_mef_session_metadata.mex*',mayo_mef)

fprintf('\n')
fprintf('Building read_mef_ts_data.mex*\n')
mex -output read_mef_ts_data ...
    read_mef_ts_data.c ...
    matmef_data.c
movefile('read_mef_ts_data.mex*',mayo_mef)

cd(mayo_mef)

% [EOF]