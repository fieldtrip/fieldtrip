function ft_omri_quality(cfg)
% function ft_omri_quality(cfg)
%
% fMRI quality assurance stack
% Make sure you have SPM8 on your path.
%
% Configuration options are
%  cfg.input            - FieldTrip buffer containing raw scans (default='buffer://localhost:1972')
%  cfg.numDummy         - how many scans to ignore initially    (default=0)
%  cfg.showRawVariation - 1 to show variation in raw scans (default), 0 to show var. in processed scans
%  cfg.clipVar          - threshold to clip variation plot with (default=300)
%  cfg.lambda           - forgetting factor for the variaton plot (default=0.9)
%  cfg.whichEcho  - which echo to process for multi-echo sequences (default = 1)
%
% (C) 2010 Stefan Klanke


%fieldtripdefs
%addpath('/home/common/matlab/spm8');

if nargin<1
	cfg = [];
end

if ~isfield(cfg, 'showRawVariation')
	cfg.showRawVariation = 1;
end

if ~isfield(cfg, 'clipVar')
	cfg.clipVar = 300;
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

if ~isfield(cfg, 'whichEcho')
	cfg.whichEcho = 1;
else
	if cfg.whichEcho < 1
		error '"whichEcho" configuration field must be >= 1';
	end
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
		warning('No protocol information found!')
		% restart loop
		pause(0.5);
		continue;
	end
	
	if cfg.whichEcho > S.numEchos
		warning('Selected echo number exceeds the number of echos in the protocol.');
		grabEcho = S.numEchos;
		fprintf(1,'Will grab echo #%i of %i\n', grabEcho, S.numEchos);
	else
		grabEcho = 1;
	end
		
	% reset motion estimates
	motEst  = [];
	maxVal  = 0;
	maxDiff = 1e-6; % zero is not possible here (for imagesc range)
	
	% Wait for numDummy scans (and drop them)
	fprintf(1,'Waiting for %i dummy samples to come in...\n', cfg.numDummy);
	while 1
		threshold = struct('nsamples', cfg.numDummy * S.numEchos);
		newNum = ft_poll_buffer(cfg.input, threshold, 500);
		if newNum.nsamples >= cfg.numDummy*S.numEchos
		   break
		end
		pause(0.01);
	end

	fprintf(1,'Starting to process\n');
	numTotal  = cfg.numDummy * S.numEchos;
	numProper = 0;
	discNum = 0;
		
	% Loop this as long as the experiment runs with the same protocol (= data keeps coming in)
	while 1
		% determine number of samples available in buffer / wait for more than numTotal
		threshold.nsamples = numTotal + S.numEchos - 1;
		newNum = ft_poll_buffer(cfg.input, threshold, 500);
		
		if newNum.nsamples < numTotal
			% scanning seems to have stopped - re-read header to continue with next trial
			break;
		end
		if newNum.nsamples < numTotal + S.numEchos
			% timeout -- go back to start of (inner) loop
			drawnow;
		    continue;
		end
		
		% this is necessary for ft_read_data
		hdr.nSamples = newNum.nsamples;
		
		index = (cfg.numDummy + numProper) * S.numEchos + grabEcho;
		fprintf('\nTrying to read %i. proper scan at sample index %d\n', numProper+1, index);
		GrabSampleT = tic;
		
		try
			% read data from buffer (only the last scan)
			dat = ft_read_data(cfg.input, 'header', hdr, 'begsample', index, 'endsample', index);
		catch
			warning('Problems reading data - going back to poll operation...');
			continue;
		end
		
		numProper = numProper + 1;
		numTotal  = numTotal + S.numEchos;
		
		maxdat = max(dat);
		if maxdat > maxVal
			maxVal = maxdat;
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
		else
			fprintf('%-30s','Registration...');
			tic; 
			[RRM, M, Mabs, procScan] = ft_omri_align_scan(RRM, rawScan); 
			toc
			motEst = [motEst; hom2six(M).*[1 1 1 180/pi 180/pi 180/pi]];
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
		
		clippedDiff = max(diffScan, cfg.clipVar);
		
		% NOTE: the following line plots "rawScan" as the current volume
		% If you'd like to see the processed scan (motion corrected) 
		% instead, exchange "rawScan" by "procScan"
		ft_omri_quality_plot(motEst, rawScan, diffScan, maxVal, cfg.clipVar);
		
		% force Matlab to update the figure
		drawnow
	end % while true	
end	
