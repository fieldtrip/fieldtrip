% Function to test a collection of MFF files

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

clear
if ispc
    baseFolder = 'z:/data/philips';
    outputFolder = '';
    error('Change folders');
elseif ismac
    baseFolder = '/Users/arno/GoogleDrive/EGI';
    outputFolder = '/data/philips/temp/';
else
    baseFolder = '/home/xxx/matlab';
end

removeICA = 0;
testtarget = 'eeglab'; % may be matlab, eeglab, fileio or fileioeeglab (compare EEGLAB and file-io) - not case sensitive

if ~strcmpi(testtarget, 'matlab') && ~strcmpi(testtarget, 'fileio')
    if ~exist('eeglab.m')
        addpath('/data/matlab/eeglab/');
    end
    eeglab; close
end

inputFilenames = { ...
    'Net Station Files from EGI/Unprocessed Continuous/32 channels/NIA_333ms_HCGSN32_test01.mff' ...
    'Net Station Files from EGI/Unprocessed Continuous/128 channels/GNG2_002_1v_cln.nsf.mff' ...
    'Net Station Files from EGI/Unprocessed Continuous/256 channels/VTD_993.1.ses.mff' ...
    'Net Station Files from EGI/Unprocessed Continuous/64 channels/NIA_333msHCGSN64_test01.mff' ...
    'Net Station Files from EGI/Continuous with Video/ZZ2_018 0416 2110.mff' ...
    'Net Station Files from EGI/Grand Average Multiple Subjects/GNG_F_Day_1_GAve.mff' ...
    'Net Station Files from EGI/Individual Averaging multiple categories/LLL_01.1_T108_0691.ave.mff' ...
    'Net Station Files from EGI/Individual Averaging multiple subjects/GNG_F_Day_1_Combined_blc_ref.mff' ...
    'Net Station Files from EGI/Segmented with multiple categories/NIA_P_013.ses.fil.seg.mff' ...
    ...
    'Net Station Files from EGI/OtherFilesArnaud/9999_20160309_011903.mff' ...
    'Net Station Files from EGI/OtherFilesArnaud/coma_LA_MOH1_20160411_060749.mff' ...
    'Net Station Files from EGI/OtherFilesArnaud/data.mff' ...
    'Net Station Files from EGI/OtherFilesArnaud/multiSubj_seg_LLL_06.cmbave.mff' ...
    'Net Station Files from EGI/OtherFilesArnaud/Nicolas_CerCo_SI_STE_test1_20151009_034300.mff' ...
    ...
    'Net Station Files from EGI/OtherFilesDavid/ATD256_1.mff' ...
    'Net Station Files from EGI/OtherFilesDavid/GNG2_014_1s_cln.seg.mff' ... % 16 file corrupted or dicontinuity to take into account (event latencies are not consistent with block latencies - this could be corrected )
    'Net Station Files from EGI/OtherFilesDavid/VTD_7Ss_bcr.gav.mff' ...
    'Net Station Files from EGI/OtherFilesDavid/VTD_993.1.ses .mff' ...
    'Net Station Files from EGI/OtherFilesDavid/4ms_5uV.nsr.mff' ...
    'Net Station Files from EGI/OtherFilesDavid/01_024 0531 1145_seg_fil_bcr_ave_WITH_AUTONOMOUS.mff' ...
    ...
    'Net Station Files from EGI/OtherFilesRobert/bug1427/Long64ChannelWithEvents.mff' ...
    'Net Station Files from EGI/OtherFilesRobert/bug1427/NS500Sine6Hz.mff' ...
    'Net Station Files from EGI/OtherFilesRobert/bug629/pilot05_test 20110120 1433.mff' ... % has autonomous data
    'Net Station Files from EGI/OtherFilesRobert/original/eeg/egi/NS500Sine6Hz.mff' ...
    'Net Station Files from EGI/Bugs/SPBI023 20150414 1357.mff' ...
    'Net Station Files from EGI/Treys_files/MMVTD_Continuous_EEG.mff' ...
    'Net Station Files from EGI/Treys_files/MMI_HC1_20180314_093330_physio_only.mff' };

comments = { ...
    '' '' '' '' '' ...
    'Differences in begintime have been verified and are OK (begintime diff less than 10 microsec)' ... % dataset 6
    '' ...
    'Differences in begintime have been verified and are OK (begintime diff less than 10 microsec)' ... % dataset 8
    '' ...
    'Differences in begintime have been verified and are OK (begintime diff less than 10 microsec)' ... % dataset 10
    '' '' '' '' '' ...
    'Differences need further investigation - has been transmitted to Philips for comments' ... % dataset 16 POTENTIAL PROBLEM
    '' '' '' ...
    'Differences in begintime have been verified and are OK (less than 10 microseconds)' ... % dataset 27
    'Differences in begintime due to version 0 vs version 3 - OK' '' ...
    'Differences in begintime due to File version 0 vs version 3 (nanoseconds vs microseconds) - OK' ...
    '' ...
    ''  ...
    'Issue with ''relativeBeginTime (event 364) has been transmitted to Philips for comments' ...
    'Differences in begintime have been verified and are OK (begintime diff than 10 microseconds)' ... % dataset 27
    };

%%
ALLEEG = [];
% the error is with dataset 7
fileDairy = sprintf('%s-%s.txt', computer, datestr(now, 'yyyymmddTHHMMSS'));
datasetsToLoad = 1:length(inputFilenames); % dataset 16 POTENTIAL PROBLEM

% scan file content
% for iFile = datasetsToLoad %1 %[20 23] %1:length(datasetsToLoad)
%     type(fullfile(fullfile(baseFolder, inputFilenames{iFile}), 'info.xml'));
%     %EEG = mff_importinfo(fullfile(baseFolder, inputFilenames{iFile}));
% end
% return

for iFile = datasetsToLoad %datasetsToLoad %1 %[20 23] %1:length(datasetsToLoad)
    errorMsg = '';
    
    fprintf('Reading file %s\n', inputFilenames{iFile});
    outputFile = fullfile(outputFolder, [ 'mffmatlabio/testexport' int2str(ispc) '_' int2str(iFile) '.mff' ] );
    if exist(outputFile)
        rmdir(outputFile, 's')
    end    
    
    if strcmpi(testtarget, 'eeglab')

        % test EEGLAB import/export
        EEG = mff_import(fullfile(baseFolder, inputFilenames{iFile}));
    
        if removeICA
            EEG = pop_eegfiltnew(EEG, [],1,826,1,[],0);
            EEG = pop_runica(EEG, 'pca', 5);
            EEG = pop_subcomp(EEG, 1);
        end
        
        mff_export(EEG, outputFile);
        EEG2 = mff_import(outputFile);
        
    else
        % test Fileio import/export
%         header = mff_fileio_read_header(fullfile(baseFolder, inputFilenames{iFile}));
%         dat    = mff_fileio_read_data(fullfile(baseFolder, inputFilenames{iFile}));
%         event  = mff_fileio_read_event(fullfile(baseFolder, inputFilenames{iFile}));
        try
            ft_defaults;        
            hdr   = ft_read_header(fullfile(baseFolder, inputFilenames{iFile}), 'headerformat', 'egi_mff_v3');
            event = ft_read_event( fullfile(baseFolder, inputFilenames{iFile}), 'eventformat', 'egi_mff_v3', 'header', hdr);
            dat   = ft_read_data( fullfile(baseFolder, inputFilenames{iFile}), 'dataformat', 'egi_mff_v3', 'header', hdr);
            EEG    = pop_fileio2(hdr, dat, event);
        catch
            l = lasterror;
            errorMsg = l.message;
        end
    
        if strcmpi(testtarget, 'fileioeeglab')
            EEG2   = mff_import(fullfile(baseFolder, inputFilenames{iFile}));
        else
            try
                ft_write_data(outputFile, dat, 'header', hdr, 'event', event, 'dataformat', 'mff');
                % mff_fileio_write(outputFile, header, dat, event);
                hdr2   = ft_read_header(outputFile, 'headerformat', 'egi_mff_v3');
                event2 = ft_read_event( outputFile, 'eventformat', 'egi_mff_v3', 'header', hdr);
                dat2   = ft_read_data(  outputFile, 'dataformat', 'egi_mff_v3', 'header', hdr);
                EEG2   = pop_fileio2(hdr2, dat2, event2);
            catch
                l = lasterror;
                errorMsg = l.message;
            end
        end
    end
    
    if ~strcmpi(testtarget, 'matlab')
        diary(fileDairy);
        disp('-------------------------')
        fprintf('File number %d\nComparing reimported file %s\n', iFile, inputFilenames{iFile});
        if isempty(errorMsg)
            eeg_compare(EEG, EEG2);
            if ~isempty(comments{iFile})
                disp(comments{iFile});
            end
        else
            disp(errorMsg);
        end
        diary
        if length(datasetsToLoad) == 1
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
            eeglab redraw; 
        end
    else
        disp('CANNOT COMPARE IMPORTED AND EXPORTED FILES USING THIS MODE');
    end
end
diary off;
return

%events = read_mff_event('../Net Station Files from EGI/Individual Averaging multiple categories/LLL_01.1_T108_0691.ave.mff', []);

%%
%EEG = mff_importsignal(inputFilenames{4});

%%

