function ft_omri_quality(cfg)

% FT_OMRI_QUALITY implements an online fMRI quality assurance stack
%
% Use as
%   ft_omri_quality(cfg)
% where cfg is a structure with configuration settings.
%
% Configuration options are
%   cfg.input            = FieldTrip buffer containing raw scans (default='buffer://localhost:1972')
%   cfg.numDummy         = how many scans to ignore initially    (default=0)
%   cfg.showRawVariation = 1 to show variation in raw scans (default), 0 to show var. in processed scans
%   cfg.clipVar          = threshold to clip variation plot with as a fraction of signal magnitude (default=0.2)
%   cfg.lambda           = forgetting factor for the variaton plot (default=0.9)
%   cfg.serial           = serial port (default = /dev/ttyS0), set [] to disable motion reporting
%   cfg.baudrate         = serial port baudrate (default = 19200)
%   cfg.maxAbs           = threshold (mm) for absolute motion before 'A' is sent to serial port, default = Inf
%   cfg.maxRel           = threshold (mm) for relative motion before 'B' is sent to serial port, default = Inf

% Copyright (C) 2010, Stefan Klanke
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ft_defaults
ft_hastoolbox('spm8', 1);

if nargin<1
  cfg = [];
end

if ~isfield(cfg, 'showRawVariation')
  cfg.showRawVariation = 1;
end

if ~isfield(cfg, 'clipVar')
  cfg.clipVar = 0.2;
end

if ~isfield(cfg, 'lambda')
  cfg.lambda = 0.9;
end

if ~isfield(cfg, 'input')
  cfg.input = 'buffer://localhost:1972';
end

if ~isfield(cfg, 'numDummy')
  cfg.numDummy = 0;
end

if ~isfield(cfg, 'baudrate')
  cfg.baudrate = 19200;
end

if ~isfield(cfg, 'serial')
  cfg.serial = '/dev/ttyS0';
end

if ~isfield(cfg, 'maxAbs')
  cfg.maxAbs = inf;
end

if ~isfield(cfg, 'maxRel')
  cfg.maxRel = inf;
end

if ~isempty(cfg.serial)
  try
    serPort = serial(cfg.serial);
    existPort = instrfind('Name',serPort.name, 'Status', 'open');
    if ~isempty(existPort)
       set(existPort, 'BaudRate', cfg.baudrate);
       serPort = existPort;
    else
       set(serPort, 'BaudRate', cfg.baudrate);
       fopen(serPort);
    end
  catch
    % the "catch me" syntax is broken on MATLAB74, this fixes it
    me = lasterror;
    serPort = [];
    disp(me.message);
  end
else
  serPort = [];
end


% Loop this forever (until user cancels)
while 1
  clear ft_read_header
  % start by reading the header from the realtime buffer
  while 1
    try
      hdr = ft_read_header(cfg.input);
      break;
    catch 
      disp('Waiting for header');
      pause(0.5);
    end
  end
  
  % Ok, we got the header, try to make sense out of it
  S = ft_omri_info_from_header(hdr);
  if isempty(S)
    ft_warning('No protocol information found!')
    % restart loop
    pause(0.5);
    continue;
  end
      
  % reset motion estimates
  motEst  = [];
  maxVal  = 0;
  maxDiff = 1e-6; % zero is not possible here (for imagesc range)
  
  % Wait for numDummy scans (and drop them)
  fprintf(1,'Waiting for %i dummy samples to come in...\n', cfg.numDummy);
  while 1
    threshold = struct('nsamples', cfg.numDummy);
    newNum = ft_poll_buffer(cfg.input, threshold, 500);
    if newNum.nsamples >= cfg.numDummy
       break
    end
    pause(0.01);
  end

  fprintf(1,'Starting to process\n');
  numTotal  = cfg.numDummy;
  numProper = 0;
  discNum = 0;
    
  % Loop this as long as the experiment runs with the same protocol (= data keeps coming in)
  while 1
    % determine number of samples available in buffer / wait for more than numTotal
    threshold.nsamples = numTotal;
    newNum = ft_poll_buffer(cfg.input, threshold, 500);
    
    if newNum.nsamples < numTotal
      % scanning seems to have stopped - re-read header to continue with next trial
      break;
    end
    if newNum.nsamples == numTotal
      % timeout -- go back to start of (inner) loop
      drawnow;
      continue;
    end
    
    % this is necessary for ft_read_data
    hdr.nSamples = newNum.nsamples;
    
    index = (cfg.numDummy + numProper) + 1;
    fprintf('\nTrying to read %i. proper scan at sample index %d\n', numProper+1, index);
    GrabSampleT = tic;
    
    try
      % read data from buffer (only the last scan)
      dat = ft_read_data(cfg.input, 'header', hdr, 'begsample', index, 'endsample', index);
    catch
      ft_warning('Problems reading data - going back to poll operation...');
      continue;
    end
    
    numProper = numProper + 1;
    numTotal  = numTotal + 1;
    
    maxdat = max(dat);
    if maxdat > maxVal
      maxVal = maxdat;
      maxDiff = cfg.clipVar * double(maxVal);
    end
    rawScan = single(reshape(dat, S.voxels));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % motion correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numProper == 1
      % set up reference model for registration
      flags = struct('mat', S.mat0);
      fprintf(1,'Setting up first non-dummy scan as reference volume...\n');
      RRM = ft_omri_align_init(rawScan, flags);
      motEst = zeros(1,6);
      procScan = single(rawScan);
                        lastPos = zeros(1,3);
    else
      fprintf('%-30s','Registration...');
      tic; 
      [RRM, M, Mabs, procScan] = ft_omri_align_scan(RRM, rawScan); 
      toc
      motPars = spm_imatrix(M);
      motEst = [motEst; motPars(1:6).*[1 1 1 180/pi 180/pi 180/pi]];
      
      if ~isempty(serPort)
        curPos = motPars(1:3);
        if any(abs(curPos) > cfg.maxAbs)
          try
            fprintf(serPort, 'A');
          catch
            % the "catch me" syntax is broken on MATLAB74, this fixes it
            me = lasterror;
            disp(me.message);
          end
          fprintf(1, 'A - too much absolute motion');
        end
        if any(abs(curPos - lastPos) > cfg.maxRel)
          try
            fprintf(serPort, 'B');
          catch
            % the "catch me" syntax is broken on MATLAB74, this fixes it
            me = lasterror;
            disp(me.message);
          end
          fprintf(1, 'B - too much relative motion');
        end
        lastPos = curPos;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % difference with previous image, after motion correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numProper == 1
      diffScan = zeros(size(procScan));
      if cfg.showRawVariation
        lastScan = rawScan;
      else
        lastScan = procScan;
      end
    else
      if cfg.showRawVariation
        % Calculation based on raw scans
        diffScan = max(cfg.lambda*diffScan, abs(rawScan - lastScan));
        lastScan = rawScan;
      else
        % Calculation based on the processed scans
        diffScan = max(cfg.lambda*diffScan, abs(procScan - lastScan));
        lastScan = procScan;
      end
      
      discNum = cfg.lambda*discNum + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % done with pre-processing, write output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fprintf('Done -- total time = %f\n', toc(GrabSampleT));
    
    clippedDiff = max(diffScan, maxDiff);
    
    % NOTE: the following line plots "rawScan" as the current volume
    % If you'd like to see the processed scan (motion corrected) 
    % instead, exchange "rawScan" by "procScan"
    ft_omri_quality_plot(motEst, rawScan, diffScan, maxVal, maxDiff);
    
    % force Matlab to update the figure
    drawnow
  end % while true  
end  
