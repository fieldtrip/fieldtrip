function image = readBESAImage(filename)

% readBESAimage reads information from exported BESA images (MSBF, BSPS,
% Sensitivity, Minimum Norm).
%
% Use as
%   image = readBESAimage(filename)
%
% The output is a structure containing the following fields:
%   Condition
%   Coordinates
%   Data
%   DataFile
%   DepthWeighting
%   Dimension
%   Frequency
%   Imagemode
%   Imagetype
%   Latency
%   Locations
%   NoiseEstimation
%   NoiseScaleFactor
%   NoiseWeighting
%   SelMeanNoise
%   Source
%   SpTmpWeighting
%   SpTmpWeightingType
%   Time
%   TimeSamples
%   Units
%   Version
%
% Modified April 26, 2006 Robert Oostenveld
% Modified November 6, 2006 Karsten Hoechstetter

if isempty(findstr(filename,'.'))
  filename = [filename,'.dat'];
end
fp = fopen(filename);

ImageVersion = fgetl(fp);

% Check Version Number
% image.Version = str2num(ImageVersion(findstr(ImageVersion,':')+1:length(ImageVersion)));

% Check type of Image
Imageinfo = '';
if ~isempty(findstr(ImageVersion,'MN'))
  image.Imagetype = 'Minimum Norm';
  image.Imagemode = '';
else
  fgetl(fp);
  ImageInfo = fgetl(fp);

  if ~isempty(findstr(ImageInfo,'Sens.'))
    image.Imagetype = 'Sensitivity';
    image.Imagemode = '';
  elseif ~isempty(findstr(ImageInfo,'MSBF'))
    image.Imagetype = 'MSBF';
    if ~isempty(findstr(ImageInfo,'Image (TF)'))
      image.Imagemode = 'Time-Frequency';
    else
      image.Imagemode = 'Time';
    end
  elseif ~isempty(findstr(ImageInfo,'MSPS'))
    image.Imagetype = 'MSPS';
    if ~isempty(findstr(ImageInfo,'Image (TF)'))
      image.Imagemode = 'Time-Frequency';
    else
      image.Imagemode = 'Time';
    end
  else
    image.Imagetype = 'Unknown';
    image.Imagemode = 'Unknown';
  end

  if ImageInfo(end-3:end) == '[dB]'
    image.Units = 'dB';
  else
    image.Units = '%';
  end
end

% Extract additional information (time and frequency, source, MN Info)
if strcmp(image.Imagemode,'Time-Frequency')
  TimeSeparator = findstr(ImageInfo,' : ');
  Blanks = findstr(ImageInfo,' ');
  [x,Index] = min(abs(Blanks-TimeSeparator));
  TimeIndex=Blanks(Index-1);
  image.Time = sscanf(ImageInfo(TimeIndex:end),'%s',4);
  image.Frequency = sscanf(ImageInfo(findstr(ImageInfo,'ms')+3:end),'%s',2);
elseif strcmp(image.Imagetype,'Sensitivity')
  image.Source = ImageInfo(findstr(ImageInfo,' - ')+3:end);
elseif strcmp(image.Imagetype,'Minimum Norm')
  fgetl(fp);
  h = fgetl(fp); image.DataFile = h(21:end);
  h = fgetl(fp); image.Condition = h(21:end);
  fgetl(fp);
  h = fgetl(fp); image.DepthWeighting = h(21:end);
  h = fgetl(fp); image.SpTmpWeighting = h(21:end);
  h = fgetl(fp); image.SpTmpWeightingType = h(21:end);
  h = fgetl(fp); image.Dimension = str2num(h(21:end));
  h = fgetl(fp); image.NoiseEstimation = h(21:end);
  h = fgetl(fp); image.NoiseWeighting = h(21:end);
  h = fgetl(fp); image.NoiseScaleFactor = str2num(h(21:end));
  h = fgetl(fp); image.SelMeanNoise = h(21:end);
  fgetl(fp);
  h = fgetl(fp); image.Locations = str2num(h(21:end));
  h = fgetl(fp); image.TimeSamples = str2num(h(21:end));
  fgetl(fp);  fgetl(fp);   fgetl(fp);
end

% Get Coordinates and Data
if ~isempty(strmatch(image.Imagetype,strvcat('MSPS','MSBF','Sensitivity')))
  % Get Coordinates
  fgetl(fp); fgetl(fp);
  h = fgetl(fp);
  hx = sscanf(h,'X: %f %f %d');
  xmin = hx(1);
  xmax = hx(2);
  xnum = hx(3);
  h = fgetl(fp);
  hy = sscanf(h,'Y: %f %f %d');
  ymin = hy(1);
  ymax = hy(2);
  ynum = hy(3);
  h = fgetl(fp);
  hz = sscanf(h,'Z: %f %f %d');
  zmin = hz(1);
  zmax = hz(2);
  znum = hz(3);
  fgetl(fp);

  image.Coordinates=struct('X',{[xmin:floor((xmax-xmin)/(xnum-1)*10000)/10000:xmax]},...
    'Y',{[ymin:floor((ymax-ymin)/(ynum-1)*10000)/10000:ymax]},...
    'Z',{[zmin:floor((zmax-zmin)/(znum-1)*10000)/10000:zmax]});

  % Get Data
  image.Data = zeros(xnum,ynum,znum);
  for z=1:znum
    fgetl(fp);
    a=fscanf(fp,'%f',[xnum,ynum]);
    for x=1:xnum
      for y=1:ynum
        image.Data(x,y,z)=a(x,y);
      end
    end
    fgetl(fp);fgetl(fp);
  end

  % Minimum Norm Image
elseif ~isempty(strmatch(image.Imagetype,('Minimum Norm')))
  h = fscanf(fp,'Latency (milliseconds):');
  image.Latency = fscanf(fp,'%f',[1,image.TimeSamples]);
  image.Coordinates = zeros(image.Locations,3);
  image.Data = zeros(image.Locations,image.TimeSamples);
  for i=1:image.Locations
    h=fscanf(fp,'%f',[1,3]);
    image.Coordinates(i,:) = h;
    image.Data(i,:)=fscanf(fp,'%f',[1,image.TimeSamples]);
  end
end

fclose(fp);
