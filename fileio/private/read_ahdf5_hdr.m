function [header] = read_ahdf5_hdr(datafile)

%read header
if nargin ~= 1
    ft_error('Wrong number of input arguments');
end

if ~isempty(datafile),
 
  id = h5readatt(datafile, '/', 'id');
  if ~strcmp(id, 'AnyWave')
      ft_error('This is not the ahf5 format.');
  end
  
  ch = h5read(datafile, '/channels');
  types = ch.type';
  labels = ch.label';
  units = ch.unit';
  refs = ch.ref';
  sr = ch.samplingRate';
  x = ch.x';
  y = ch.y';
  z = ch.z';
  ox = ch.ox';
  oy = ch.oy';
  oz = ch.oz';
  
  header.numberOfBlocks = h5readatt(datafile, '/Blocks', 'numberOfBlocks');
  header.numberOfSamples = h5readatt(datafile, '/Blocks/1', 'numberOfSamples');
  for i=1:size(types, 1)
      header.channels(i).samplingRate = sr(i);
      header.channels(i).type = strcat(types(i, 1:16));
      switch header.channels(i).type
          case 'SEEG'
              header.type{i} = 'seeg';
          case 'EEG'
              header.type{i} = 'eeg';
          case 'MEG'
              header.type{i} = 'meg';
          case 'ECG'
            header.type{i} = 'ecg';
          case 'EMG'
            header.type{i} = 'emg';
      end
                 
      header.channels(i).label = strcat(labels(i, 1:16));
      header.channels(i).ref = strcat(refs(i, 1:16));
      header.channels(i).unit = strcat(units(i, 1:4));
      header.channels(i).x = x(i);
      header.channels(i).y = y(i);
      header.channels(i).z = z(i);
      header.channels(i).ox = ox(i);
      header.channels(i).oy = oy(i);
      header.channels(i).oz = oz(i);     
      % create label cell-array
      header.label{i} =  header.channels(i).label;
      % create reference cell-array
      header.reference{i} = header.channels(i).ref;
      % create chantype
      if  strcmp(header.channels(i).unit, '�V')
          header.unit{i} = 'uV';
      else
          header.unit{i} = header.channels(i).unit;
      end
  end 
end
