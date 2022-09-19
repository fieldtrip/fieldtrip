function test_nicolet_reading_onefile(path_to_load,nicoletFileName,asciiVerifyFileName,Fs,nChans,nSamples,nTrials,startDateTime)

datetimeTolerance = 0.00000000001;
asciiDataTolerance = 0.01;
fullNicoletFile1 = fullfile(path_to_load,nicoletFileName);
fullasciiVerifyFile = fullfile(path_to_load,asciiVerifyFileName);

disp(' ');
disp('Testing Nicolet/Nervus data import with file:');
disp(fullNicoletFile1);
disp('=====');

try
  disp('Reading data file header...');  
  hdr = ft_read_header(fullNicoletFile1);
catch ME
  disp(ME);
  error('Failed reading header of Nicolet/Nervus EEG file');
end

try
  disp('Checking header');
  assert(hdr.Fs==Fs,['Sampling rate should be ' num2str(Fs) ', was ' num2str(hdr.Fs)]);
  assert(hdr.nChans==nChans,['Channel count should be ' num2str(nChans) ', was ' num2str(hdr.nChans)]);
  assert(hdr.nSamples==nSamples,['Sample count should be ' num2str(nSamples) ', was ' num2str(hdr.nSamples)]);
  assert(hdr.nTrials==nTrials,['Trial count should be ' num2str(nTrials) ', was ' num2str(hdr.nTrials)]);
  assert(hdr.orig.startDateTime-startDateTime<datetimeTolerance,['startDateTime should be ' num2str(datestr(startDateTime)) ',was ' num2str(datestr(hdr.orig.startDateTime))]);
  assert(strcmp(hdr.orig.reference,'common'),'Reference should be common');
catch ME
  disp(ME);
  error(['Failed verifying header of Nicolet/Nervus EEG: ' ME.message]);
end

try
  disp('Reading data ...');
  dataopts = {};
  alldata = ft_read_data(fullNicoletFile1, 'header', hdr, dataopts{:});
catch ME
  disp(ME);
  error('Failed reading data from Nicolet/Nervus EEG file');
end
try
  assert(size(alldata,1)==nChans,['Channel data should be ' num2str(nChans) ' long']);
  assert(size(alldata,2)==nSamples,['Sample data should be ' num2str(nSamples) ' long']);    
catch ME
  disp(ME);
  error('Failed verifying data size from Nicolet/Nervus EEG file');
end

if (isempty(asciiVerifyFileName))
  warning('Skipping ASCII export verification')
else
  try 
    disp('Reading ASCII-exported data to verify signals');
    fid = fopen(fullasciiVerifyFile,'rt');
    textHeader = strtrim(regexp(fgetl(fid),'\t','split'));  
    fclose(fid);
    textData = dlmread(fullasciiVerifyFile,'\t',1,0)';

  catch ME
    disp(ME);
    error('Failed reading verification ASCII data file');
  end
  
  try 
  disp('Comparing ASCII-exported data to Nicolet data');
  assert(size(alldata,1)<=size(textData,1),'Nicolet data should have same or fewer columns than text data');
  assert((size(alldata,2))==size(textData,2),'Nicolet data should have same sample count as text data');

  for col=1:size(alldata,1)
    for sample=1:(size(alldata,2)-1)
      assert(abs(alldata(col,sample)-textData(col,sample))<asciiDataTolerance,'Sample inequality found');
    end
  end
  disp('Data are equal up to the limits of Nicolet ASCII export');

  catch ME
    disp(ME);
    error('Failed comparing text data to Nicolet data');
  end
end
  
end