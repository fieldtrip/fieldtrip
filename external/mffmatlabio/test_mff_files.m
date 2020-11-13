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

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    confirm_recursive_rmdir(false);
end

removeICA = 0;
testtarget = 'eeglab'; % may be matlab, eeglab, fileio or fileioeeglab (compare EEGLAB and file-io) - not case sensitive

if ~strcmpi(testtarget, 'matlab') && ~strcmpi(testtarget, 'fileio')
    if ~exist('eeglab.m')
        addpath('/data/matlab/eeglab/');
    end
    eeglab; close
end

inputFilenames = [];
inputFilenames(end+1).file = 'MFF_Files/Bugs/RHYM_MFF_issue19/AA_1756.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/RHYM_MFF_issue19/AA_1845.mff';

inputFilenames(end+1).file = 'MFF_Files/Unprocessed Continuous/32 channels/NIA_333ms_HCGSN32_test01.mff';
inputFilenames(end+1).file = 'MFF_Files/Unprocessed Continuous/128 channels/GNG2_002_1v_cln.nsf.mff';
inputFilenames(end+1).file = 'MFF_Files/Unprocessed Continuous/256 channels/VTD_993.1.ses.mff';
inputFilenames(end+1).file = 'MFF_Files/Unprocessed Continuous/64 channels/NIA_333msHCGSN64_test01.mff';
inputFilenames(end+1).file = 'MFF_Files/Continuous with Video/ZZ2_018 0416 2110.mff';
inputFilenames(end+1).file = 'MFF_Files/Grand Average Multiple Subjects/GNG_F_Day_1_GAve.mff'; inputFilenames(end).info = 'Differences in begintime have been verified and are OK (begintime diff less than 10 microsec)'; % dataset 6
inputFilenames(end+1).file = 'MFF_Files/Individual Averaging multiple categories/LLL_01.1_T108_0691.ave.mff';
inputFilenames(end+1).file = 'MFF_Files/Individual Averaging multiple subjects/GNG_F_Day_1_Combined_blc_ref.mff'; inputFilenames(end).info =  'Differences in begintime have been verified and are OK (begintime diff less than 10 microsec)'; % dataset 8
inputFilenames(end+1).file = 'MFF_Files/Segmented with multiple categories/NIA_P_013.ses.fil.seg.mff';

inputFilenames(end+1).file = 'MFF_Files/OtherFilesArnaud/9999_20160309_011903.mff'; inputFilenames(end).info = 'Differences in begintime have been verified and are OK (begintime diff less than 10 microsec)'; % dataset 10
inputFilenames(end+1).file = 'MFF_Files/OtherFilesArnaud/coma_LA_MOH1_20160411_060749.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesArnaud/data.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesArnaud/multiSubj_seg_LLL_06.cmbave.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesArnaud/Nicolas_CerCo_SI_STE_test1_20151009_034300.mff';

inputFilenames(end+1).file = 'MFF_Files/OtherFilesDavid/ATD256_1.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesDavid/GNG2_014_1s_cln.seg.mff'; inputFilenames(end).info = 'Differences need further investigation - has been transmitted to Philips for comments'; % dataset 16 POTENTIAL PROBLEM - % 16 file corrupted or dicontinuity to take into account (event latencies are not consistent with block latencies - this could be corrected )
inputFilenames(end+1).file = 'MFF_Files/OtherFilesDavid/VTD_7Ss_bcr.gav.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesDavid/VTD_993.1.ses .mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesDavid/4ms_5uV.nsr.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesDavid/01_024 0531 1145_seg_fil_bcr_ave_WITH_AUTONOMOUS.mff';

inputFilenames(end+1).file = 'MFF_Files/OtherFilesRobert/bug1427/Long64ChannelWithEvents.mff'; inputFilenames(end).info = 'Differences in begintime have been verified and are OK (less than 10 microseconds)';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesRobert/bug1427/NS500Sine6Hz.mff'; inputFilenames(end).info = 'Differences in begintime due to version 0 vs version 3 - OK';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesRobert/bug629/pilot05_test 20110120 1433.mff'; % has autonomous data
inputFilenames(end+1).file = 'MFF_Files/OtherFilesRobert/original/eeg/egi/NS500Sine6Hz.mff';  'Differences in begintime due to File version 0 vs version 3 (nanoseconds vs microseconds) - OK';

inputFilenames(end+1).file = 'MFF_Files/Bugs/SPBI023_20150414_1357.mff';
inputFilenames(end+1).file = 'MFF_Files/Treys_files/MMVTD_Continuous_EEG.mff';
inputFilenames(end+1).file = 'MFF_Files/Treys_files/MMI_HC1_20180314_093330_physio_only.mff'; inputFilenames(end).info = 'Issue with ''relativeBeginTime (event 364 differs) - field is not exported back and is not 0 in original file - field not used anyway';

inputFilenames(end+1).file = 'MFF_Files/OtherFilesJoe/SLI_30.ave.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesJoe/EEG-fMRI_seg_1.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesJoe/1001_fil_bcr_ref_seg.mff';

inputFilenames(end+1).file = 'MFF_Files/OtherFilesKyle/Grand_Average_Example.mff'; inputFilenames(end).info = 'Differences in MFF keys and begin latencies have been verified - latencies are identical by the way';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesKyle/Individual_Subject_Segmented_Example.mff'; inputFilenames(end).info = 'Latency differences have been verified';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesKyle/Multiple_Subject_Averages_Example.mff';

inputFilenames(end+1).file = 'MFF_Files/OtherFilesMark/MMVTD_fil_fil_seg_1_ave.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesMark/MMVTD_fil_fil_seg_1.mff';
inputFilenames(end+1).file = 'MFF_Files/OtherFilesMark/MMVTD_fil.mff';

inputFilenames(end+1).file = 'MFF_Files/Bugs/14_20191101_105925.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/catherine/PS38_V1_S1_fMRI_Vis_block_20180516_104228.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/catherine/TestIO_PS38_V1_S1_fMRI_Vis_block_20180516_104228_mrr_fil_bcr_qrs_obs_seg.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/catherine/TestIO_PS38_V1_S1_fMRI_Vis_block_20180516_104228.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/forArnoFromJoe2/demo_20130911_1214NS5_seg.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/forArnoFromJoe2/GNG2_002_1v_cln.seg.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/forArnoFromJoe2/NIA_P_013.ses.fil.seg.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue1/samp_file.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue4/0027_Sungjin_20170403_114721_Hfil_Lfil_seg_bcr_blc_ref.mff'; inputFilenames(end).info = 'File with categories and also missing data portions - cannot be properly handled by EEGLAB';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue6/Baseline_PreNback.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue16/s34_control_20190507_022838.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue17/AV1_AV.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/KHA_SEP_20180529_110628_epoched_300ms_ave.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/KHA_SEP_20180529_110628_epoched_300ms.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/KHA_SEP_20180529_110628.mff'; inputFilenames(end).info = 'One channel difference OK (E1001 vs E129) because eeg_importcoordinates does not store chan number (1001 in original)';
inputFilenames(end+1).file = 'MFF_Files/Bugs/mattan/2002 20120824 dis copy.mff'; inputFilenames(end).info = 'Differences in begintime have been verified and are OK (begintime diff than 10 microseconds)';
inputFilenames(end+1).file = 'MFF_Files/Bugs/Piezo_EqptTest1K_8_6_15.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/Priefer_AV1_AV_import_export_issue_to_test_in_NETstation/AV1_AV.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/SPBI023_20150414_1357.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/trevor/SSDE_P004_RSEEG_EO_BS_20180922_042706.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/zhaodi.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue21/EEG_VisualMotor_Run01.mff';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue22/CNTL001_P50Sup_T1_20150608_1530.mff'; inputFilenames(end).info = 'Differences checked and OK - reference renamed for some reason';
inputFilenames(end+1).file = 'MFF_Files/Bugs/issue25/2197LA_creativity1.mff'; inputFilenames(end).info = 'Differences checked and OK - reference renamed for some reason';

%%
ALLEEG = [];
% the error is with dataset 7
fileDairy = sprintf('%s-%s.txt', computer, datestr(now, 'yyyymmddTHHMMSS'));
datasetsToLoad = 1:length(inputFilenames); % dataset 16 POTENTIAL PROBLEM

% scan file content
% for iFile = datasetsToLoad %1 %[20 23] %1:length(datasetsToLoad)
%     type(fullfile(fullfile(baseFolder, inputFilenames(iFile).file), 'info.xml'));
%     %EEG = mff_importinfo(fullfile(baseFolder, inputFilenames{iFile}));
% end
% return

for iFile = length(datasetsToLoad)
    errorMsg = '';
    
    fprintf('Reading file %s\n', inputFilenames(iFile).file);
    outputFile = fullfile(outputFolder, [ 'mffmatlabio/testexport' int2str(ispc) '_' int2str(iFile) '.mff' ] );
    if exist(outputFile)
        rmdir(outputFile, 's')
    end    
    
    if strcmpi(testtarget, 'eeglab')

        % test EEGLAB import/export
        EEG = mff_import(fullfile(baseFolder, inputFilenames(iFile).file));
    
        if removeICA
            EEG = pop_eegfiltnew(EEG, [],1,826,1,[],0);
            EEG = pop_runica(EEG, 'pca', 5);
            EEG = pop_subcomp(EEG, 1);
        end
        
        mff_export(EEG, outputFile);
        EEG2 = mff_import(outputFile);
        
    else
        
        if strcmpi(testtarget, 'matlab')
            EEG2   = mff_import(fullfile(baseFolder, inputFilenames(iFile).file));
        else
            
            try
                ft_defaults;
                hdr   = ft_read_header(fullfile(baseFolder, inputFilenames(iFile).file), 'headerformat', 'egi_mff_v3');
                event = ft_read_event( fullfile(baseFolder, inputFilenames(iFile).file), 'eventformat', 'egi_mff_v3', 'header', hdr);
                dat   = ft_read_data(  fullfile(baseFolder, inputFilenames(iFile).file), 'dataformat', 'egi_mff_v3', 'header', hdr);
                EEG    = pop_fileio2(hdr, dat, event);
            catch
                l = lasterror;
                errorMsg = l.message;
            end
            
            if strcmpi(testtarget, 'fileioeeglab')
                EEG2   = mff_import(fullfile(baseFolder, inputFilenames(iFile).file));
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
    end
    
    if ~strcmpi(testtarget, 'matlab')
        diary(fileDairy);
        disp('-------------------------')
        fprintf('File number %d\nComparing reimported file %s\n', iFile, inputFilenames(iFile).file);
        if isempty(errorMsg)
            eeg_compare(EEG, EEG2);
            if ~isempty(inputFilenames(iFile).info)
                disp(inputFilenames(iFile).info);
                disp(inputFilenames(iFile).info);
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

%events = read_mff_event('../MFF_Files/Individual Averaging multiple categories/LLL_01.1_T108_0691.ave.mff', []);

%%
%EEG = mff_importsignal(inputFilenames{4});

%%

