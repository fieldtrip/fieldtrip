function test_nicolet_reading

% MEM 2gb
% WALLTIME 00:30:00
% DEPENDENCY test_nicolet_reading_onefile read_nervus_header

% function to test reading of Nicolet/Nervus EEG files
% one is a file from 2006 (older Nicolet format)
% one is a file from 2018 (newer Nicolet format)

%path_to_load = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet');
path_to_load = 'C:\Midlertidig_Lagring\nicolet-test-data';

%%
file1 = '2018_Patient59_EEG-OPPTAKER-1_t1-NICOLET.e';
file1ascii = '2018_Patient59_EEG-OPPTAKER-1_t1-ASCII.txt';
test_nicolet_reading_onefile(path_to_load,file1,file1ascii,500,24,45501,1,datetime(2018,06,18,08,07,24));



%%
file2 = '2006_janbrogger.e';
file2ascii = '2006_janbrogger.txt';
test_nicolet_reading_onefile(path_to_load,file2,file2ascii,256,24,307201,1,datetime(2006,06,09,13,42,41));

%%
%file from 2020
file3 = '2020_Patient4_EEG246_t1.e';
file3ascii = '2020_Patient4_EEG246_t1.txt';
test_nicolet_reading_onefile(path_to_load,file3,file3ascii,500,27,605501,1,datetime(2020,03,18,14,20,34));

%%
%file from 2012
file4 = '2012_Patient2_EEG246_t1.e';
file4ascii = '';
test_nicolet_reading_onefile(path_to_load,file4,file4ascii,500,28,714001,1,datetime(2018,06,18,08,07,24));

%%
% Load a list of 10 random EEGs per year from 1997 to 2020, and try to
% read:
% SELECT YEAR(CAST(CASE WHEN dbo.tblTest.dodtRecordingStartTime > 0 THEN dbo.tblTest.dodtRecordingStartTime-2.0  ELSE  2*CAST(dbo.tblTest.dodtRecordingStartTime AS INT) - 2.0 +  ABS(dbo.tblTest.dodtRecordingStartTime) END as datetime)) AS Year,
% 	dbo.eegfiles.strName AS FileName, 
% 	dbo.eegfiles.strPath AS FilePath
% FROM    dbo.tblStudy 
%     INNER JOIN dbo.tblStudyTest ON dbo.tblStudy.guidStudyID = dbo.tblStudyTest.guidStudyID 
%     INNER JOIN dbo.tblTest ON dbo.tblStudyTest.guidTestID = dbo.tblTest.guidTestID 
%     INNER JOIN dbo.eegfiles ON dbo.tblTest.lTestID = dbo.eegfiles.lTest_Id     
% ORDER BY dbo.tblTest.dodtRecordingStartTime, NEWID()
testlist = importdata('C:\Midlertidig_Lagring\nicolet-test-data\randomsamplebyyear.txt','\t',1); 
testresult = struct();
testresult.year = 0;
testresult.testcount = 0;
testresult.successcount = 0;
startYear = 2000;
for i=2:size(testlist.textdata,1)
    year = str2double(testlist.textdata(i,1));
        
    testresult(year-startYear ).year = year;
    if isempty(testresult(year-startYear ).testcount)
        testresult(year-startYear ).testcount = 0;
    end
    testresult(year-startYear ).testcount = testresult(year-startYear ).testcount+1;
    
    fileName = char(testlist.textdata(i,4));
    filePath = char(testlist.textdata(i,5));    
    fullfile1 = fullfile(filePath,fileName);    
    %disp(fullfile1);
    success = 0;
    try 
        if isempty(testresult(year-startYear ).successcount)
            testresult(year-startYear ).successcount = 0;
        end
        
        hdr = ft_read_header(fullfile1);
        dataopts = {};
        %ft_read_data(fullfile1, 'header', hdr, dataopts{:});
        success = 1;        
        testresult(year-startYear ).successcount = testresult(year-startYear ).successcount+1;
    catch ME
        disp(ME);
    end
    disp(['Read file number ' num2str(i-1) ' from year ' num2str(year) ' success: ' num2str(success)]);        
end
struct2table(testresult)
    

%%
% Test code requested by Robert. For file 1
disp('Checking that files can be displayed. File 1');
filepath1 = fullfile(path_to_load,file1);

hdr = ft_read_header(filepath1);
event = ft_read_event(filepath1);
dataopts = {};
alldata = ft_read_data(filepath1, 'header', hdr, dataopts{:});

% open databrowser for file 1
cfg            = [];
cfg.dataset    = filepath1;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);

%%
% Tests code requested by Robert. For file 2
filepath2 = fullfile(path_to_load,file2);

hdr = ft_read_header(filepath2);
event = ft_read_event(filepath2);
dataopts = {};
alldata = ft_read_data(filepath2, 'header', hdr, dataopts{:});

% open databrowser for file 2
cfg            = [];
cfg.dataset    = filepath2;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);
