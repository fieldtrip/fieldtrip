% [data, hdr] = opensvl(filename)
%
% Reads a CTF SAM (.svl) file.

function [data, hdr] = read_ctf_svl(filename)
  
  fid = fopen(filename, 'rb', 'ieee-be', 'ISO-8859-1');
  
  if fid <= 0
    error('Could not open SAM file: %s\n', filename);
  end
  

  % ----------------------------------------------------------------------
  % Read header.
  hdr.identity = fread(fid, 8, '*char')'; % 'SAMIMAGE'
  hdr.version = fread(fid, 1, 'int32'); % SAM file version.
  hdr.setName = fread(fid, 256, '*char')'; % Dataset name.
  hdr.numChans = fread(fid, 1, 'int32');
  hdr.numWeights = fread(fid, 1, 'int32'); % 0 for static image.
  if(hdr.numWeights ~= 0)
    warning('hdr.numWeights ~= 0');
  end
  
  fread(fid,1,'int32'); % Padding to next 8 byte boundary.
  
  hdr.xmin = fread(fid, 1, 'double'); % Bounding box coordinates (m).
  hdr.xmax = fread(fid, 1, 'double');
  hdr.ymin = fread(fid, 1, 'double');
  hdr.ymax = fread(fid, 1, 'double');
  hdr.zmin = fread(fid, 1, 'double');
  hdr.zmax = fread(fid, 1, 'double');
  hdr.stepSize = fread(fid, 1, 'double'); % m
  
  hdr.hpFreq = fread(fid, 1, 'double'); % High pass filtering frequency (Hz).
  hdr.lpFreq = fread(fid, 1, 'double'); % Low pass.
  hdr.bwFreq = fread(fid, 1, 'double'); % Bandwidth
  hdr.meanNoise = fread(fid, 1, 'double'); % Sensor noise (T).
  
  hdr.mriName = fread(fid, 256, '*char')';
  hdr.fiducial.mri.nas = fread(fid, 3, 'int32'); % CTF MRI voxel coordinates?
  hdr.fiducial.mri.rpa = fread(fid, 3, 'int32');
  hdr.fiducial.mri.lpa = fread(fid, 3, 'int32');
  
  hdr.SAMType = fread(fid, 1, 'int32'); % 0: image, 1: weights array, 2: weights list.
  hdr.SAMUnit = fread(fid, 1, 'int32'); 
  % Possible values: 0 coefficients Am/T, 1 moment Am, 2 power (Am)^2, 3 Z,
  % 4 F, 5 T, 6 probability, 7 MUSIC.
  
  fread(fid, 1, 'int32'); % Padding to next 8 byte boundary.
  
  if hdr.version > 1
    % Version 2 has extra fields.
    hdr.fiducial.head.nas = fread(fid, 3, 'double'); % CTF head coordinates?
    hdr.fiducial.head.rpa = fread(fid, 3, 'double');
    hdr.fiducial.head.lpa = fread(fid, 3, 'double');
    hdr.SAMUnitName = fread(fid, 32, '*char')';
  % Possible values: 'Am/T' SAM coefficients, 'Am' source strength,
  % '(Am)^2' source power, ('Z', 'F', 'T') statistics, 'P' probability.
  end
  
  
  % ----------------------------------------------------------------------
  % Read image data.
  data = fread(fid, inf, 'double'); 
  fclose(fid);
  
  % Raw image data is ordered as a C array with indices: [x][y][z], meaning
  % z changes fastest and x slowest.  These x, y, z axes point to ALS
  % (anterior, left, superior) respectively in real world coordinates,
  % which means the voxels are in SLA order.

  
  % ----------------------------------------------------------------------
  % Post processing.
  
  % Change from m to mm.
  hdr.xmin = hdr.xmin * 1000;
  hdr.ymin = hdr.ymin * 1000;
  hdr.zmin = hdr.zmin * 1000;
  hdr.xmax = hdr.xmax * 1000;
  hdr.ymax = hdr.ymax * 1000;
  hdr.zmax = hdr.zmax * 1000;
  hdr.stepSize = hdr.stepSize * 1000;

  % Number of voxels in each dimension.
  hdr.dim = [round((hdr.xmax - hdr.xmin)/hdr.stepSize) + 1, ...
    round((hdr.ymax - hdr.ymin)/hdr.stepSize) + 1, ...
    round((hdr.zmax - hdr.zmin)/hdr.stepSize) + 1];
  
  data = reshape(data, hdr.dim([3, 2, 1]));
  
  % Build transformation matrix from raw voxel coordinates (indexed from 1)
  % to head coordinates in mm.  Note that the bounding box is given in
  % these coordinates (in m, but converted above).
  % Apply scaling.
  hdr.transform = diag([hdr.stepSize * ones(1, 3), 1]);
  % Reorder directions.
  hdr.transform = hdr.transform(:, [3, 2, 1, 4]);
  % Apply translation.
  hdr.transform(1:3, 4) = [hdr.xmin; hdr.ymin; hdr.zmin] - hdr.stepSize;
  % -step is needed since voxels are indexed from 1.
    
end

