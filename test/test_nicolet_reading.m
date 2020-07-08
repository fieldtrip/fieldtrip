function test_nicolet_reading

% MEM 2gb
% WALLTIME 00:30:00
% DEPENDENCY test_nicolet_reading_onefile read_nervus_header

% function to test reading of Nicolet/Nervus EEG files
% one is a file from 2006 (older Nicolet format)
% one is a file from 2018 (newer Nicolet format)

path_to_load = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet');
%path_to_load = 'C:\Midlertidig_Lagring\FieldTripNicoletTestData';

%%
% file1 = '2018_Patient59_EEG-OPPTAKER-1_t1-NICOLET.e';
% file1ascii = '2018_Patient59_EEG-OPPTAKER-1_t1-ASCII.txt';
% %file1 = 'Patient59_EEG-OPPTAKER-1_t1-NICOLET.e';
% %file1ascii = 'Patient59_EEG-OPPTAKER-1_t1-ASCII.txt';
% test_nicolet_reading_onefile(path_to_load,file1,file1ascii,500,24,45500,1,datetime(2018,06,18,08,07,24));

%%
file2 = '2006_janbrogger.e';
file2ascii = '2006_janbrogger.txt';
%file2 = 'janbrogger.e';
%file2ascii = 'janbrogger.txt';
test_nicolet_reading_onefile(path_to_load,file2,file2ascii,256,24,307200,1,datetime(2006,06,09,13,42,41));

%%
% Test code requested by Robert. For file 1
% disp('Checking that files can be displayed. File 1');
% filepath1 = fullfile(path_to_load,file1);
% 
% hdr = ft_read_header(filepath1);
% event = ft_read_event(filepath1);
% dataopts = {};
% alldata = ft_read_data(filepath1, 'header', hdr, dataopts{:});

% open databrowser for file 1
cfg            = [];
%cfg.dataset    = filepath1;
cfg.continuous = 'yes';
cfg.channel    = 'all';
%data           = ft_preprocessing(cfg);
cfg.viewmode   = 'vertical';
%ft_databrowser(cfg, data);

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


%%%%%%%%%%%%%
% Tests Nicolet reading on clinical EEGs - not available for regular
% fieldtrip testingsting
%%%%%%%%%%%%%
%
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

%%%%%%%%%%%%%%
% MATLAB test code that uses the SQL query results above
%%%%%%%%%%%%%%
% testlist = importdata('C:\Midlertidig_Lagring\FieldTripNicoletTestData\randomsamplebyyear.txt','\t',1); 
% testresult = struct();
% testresult.year = 0;
% testresult.testcount = 0;
% testresult.successcount = 0;
% startYear = 1997;
% for i=2:size(testlist.textdata,1)
%     tic
%     year = str2double(testlist.textdata(i,1));
%         
%     testresult(year-startYear ).year = year;
%     if isempty(testresult(year-startYear ).testcount)
%         testresult(year-startYear ).testcount = 0;
%     end
%     testresult(year-startYear ).testcount = testresult(year-startYear ).testcount+1;
%     
%     fileName = char(testlist.textdata(i,4));
%     filePath = char(testlist.textdata(i,5));    
%     fullfile1 = fullfile(filePath,fileName);    
%     %disp(fullfile1);
%     success = 0;
%     try 
%         if isempty(testresult(year-startYear ).successcount)
%             testresult(year-startYear ).successcount = 0;
%         end
%         
%         hdr = ft_read_header(fullfile1);
%         dataopts = {};
%         ft_read_data(fullfile1, 'header', hdr, dataopts{:});
%         success = 1;        
%         testresult(year-startYear ).successcount = testresult(year-startYear ).successcount+1;
%     catch ME
%         disp(ME);
%     end
%     disp(['Read file number ' num2str(i-1) ' from year ' num2str(year) ' success: ' num2str(success)]);        
%     toc
% end
% struct2table(testresult)


% % 
%  Test result from 10 clinical EEGs per calendar year 1998-2020
%     year    testcount    successcount
%     ____    _________    ____________
% 
%     1998       10              0     
%     1999       10              0     
%     2000       10              0     
%     2001       10              0     
%     2002       10              0     
%     2003       10              0     
%     2004       10              0     
%     2005       10              0     
%     2006       10              0     
%     2007       10              0     
%     2008       10              0     
%     2009       10              0     
%     2010       10              0     
%     2011       10              0     
%     2012       10             10     
%     2013       10             10     
%     2014       10             10     
%     2015       10             10     
%     2016       10             10     
%     2017       10             10     
%     2018       10             10     
%     2019       10             10     
%     2020       10             10  
    

