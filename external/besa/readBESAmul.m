function mul = readBESAmul(filename)

% readBESAmul read information from a *.mul file
%
% Use as
%   mul = readBESAMul(filename)
%
% The output is a structure containing the following fields:
%   Npts: number of sample ponts
%   TSB: latency of the first sample in the data files
%   DI: time interval between two sample points
%   Scale: scaling factor
%   ChannelLabels: Channel labels
%   Data: data matrix [Npts x Number of Channels]
%
% Last modified April 26, 2006 Robert Oostenveld

if isempty(findstr(filename,'.'))
  filename = [filename,'.mul'];
end

fp = fopen(filename);

% Read the file
if (fp)
  % Read header of .mul file
  headline=fscanf(fp,'TimePoints=%f  Channels=%f  BeginSweep[ms]=%f  SamplingInterval[ms]=%f  Bins/uV=%f');
  mul.Npts=headline(1);
  NChan=headline(2);
  mul.TSB=headline(3);
  mul.DI=headline(4);
  mul.Scale=headline(5);
  fgets(fp);
  for Channel=1:NChan
    mul.ChannelLabels(Channel) = cellstr(fscanf(fp,'%s ',1));
  end

  mul.data=zeros(mul.Npts,NChan);
  for i=1:mul.Npts
    mul.data(i,:)=fscanf(fp,'%f',[1 NChan]);
  end
end
fclose(fp);
