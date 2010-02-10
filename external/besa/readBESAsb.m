function [time,data,nEpochs] = readBESAsb(filename)

% readBESAsb reads information from a *.dat, i.e. a simple binary data file
%
% This function requires a file with the same basename as the data name but
% suffix '.generic' or '.gen' to be located in the same folder. This file is 
% generated automatically by BESA during file export.
%
% Use as
%   [time,data,nEpochs] = readBESAsb(filename)
%
% The following output is generated:
%   time: a vector containing the time points corresponding to the data
%         points. Negative times are prestimulus
%   data: The data matrix with dimension [nChannels x nEpochs x nSamples],
%         where nChannels is the number of Channels and nSamples the number
%         of samples within one epoch
%   nEpochs (optional): The number of epochs contained in the file

% Modified February 21, 2007 Karsten Hoechstetter
% Modified April 24, 2007 Robert Oostenveld
% Modified September 24, 2009 Karsten Hoechstetter

if isempty(findstr(filename,'.dat'))
  filename = [filename,'.dat'];
end

fid=fopen([filename(1:end-4),'.generic'],'r');
if fid==-1
  fid=fopen([filename(1:end-4),'.gen'],'r');
end

fscanf(fid,'BESA Generic Data\n');
nChannels = fscanf(fid,'nChannels=%i\n');
sRate = fscanf(fid,'sRate=%f\n');
nSamples = fscanf(fid,'nSamples=%i\n');
format = fscanf(fid,'format=%s');
file = fscanf(fid,'\nfile=%s');
prestimulus = fscanf(fid,'prestimulus=%f\n');   if isempty(prestimulus),prestimulus=0;end
epochs = fscanf(fid,'epochs=%i\n');             if isempty(epochs),epochs=1;end
fclose(fid);

time=[-prestimulus:1/sRate*1000:(nSamples/epochs-1)*1000/sRate-prestimulus];

fid=fopen(filename, 'r', 'ieee-le');
xdata=fread(fid,[nChannels,nSamples],'float32');
fclose(fid);

data=zeros(nChannels,epochs,nSamples/epochs);
for i=1:epochs
  data(:,i,:)=xdata(:,1+(i-1)*nSamples/epochs:i*nSamples/epochs);
end
nEpochs=epochs;

