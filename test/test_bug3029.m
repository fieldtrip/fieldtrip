function test_bug3029

% MEM 5gb
% WALLTIME 00:20:00

% TEST test_bug3029 ft_sourceanalysis prepare_freq_matrices

% this test function tests the correctness of the data handling in
% ft_sourceanalysis, with respect to:
% -refchan handling when using DICS (with refchan)
% -csd-matrix handling when using various versions of keeptrials in cfg

filename = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg'), 'freq_mtmfft_fourier_trl_ctf275');
load(filename);

filename = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/vol'), 'Subject01vol_singleshell');
load(filename);
vol = ft_datatype_headmodel(vol);
[vol.bnd.nrm] = normals(vol.bnd.pos,vol.bnd.tri);

cfg           = [];
cfg.headmodel = vol;
cfg.grad      = freq.grad;
cfg.grid.resolution = 2.5;
s             = ft_prepare_sourcemodel(cfg);

cfg           = [];
cfg.headmodel = vol;
cfg.grad      = freq.grad;
cfg.grid      = s;
cfg.channel   = 'MEG';
lf            = ft_prepare_leadfield(cfg);

% With 3029 unfixed the following cfg gives a reshape problem when averaging
% the Cf, where the assumption was that it was an NtrialxNchanxNchan matrix,
% whereas Cf already was NchanxNchan
if 0,
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
s1            = ft_sourceanalysis(cfg, freq);
end

% the following runs through irrespective of the bug being fixed, yet the
% behavior might be unexpected: no single trial output is generated.
if 0,
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.keeptrials = 'yes';
s2             = ft_sourceanalysis(cfg, freq);
end

if 0,
% with the bug unfixed, this gives a problem with the refchan handling, due
% to ft_selectdata only keeping the MEG channels from an early stage
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
s3             = ft_sourceanalysis(cfg, freq);

cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
cfg.keeptrials = 'yes';
s4             = ft_sourceanalysis(cfg, freq);

assert(isalmostequal(s3.avg.coh, s4.avg.coh, 'reltol', 1e-8));
end

% Tthe following cfg should give an error due to the mismatch in
% channels in data and leadfield
if 0,
try,
	cfg           = [];
	cfg.method    = 'dics';
	cfg.frequency = 10;
	cfg.headmodel = vol;
	cfg.grid      = lf;
	s5            = ft_sourceanalysis(cfg, freq);
catch
	% this should be solved with proper error handling
	if ~isempty(strfind(lasterr,'There''s a mismatch'))
		% this is OK
	else
		error('there''s a problem');
	end
end
end

if 0,
% check whether it also works with pcc as a method
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
cfg.dics.keepcsd = 'yes';
cfg.dics.realfilter = 'yes';
s3             = ft_sourceanalysis(cfg, freq);

cfg           = [];
cfg.method    = 'pcc';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
s3b           = ft_sourceanalysis(cfg, freq);
s3c           = ft_sourcedescriptives([], s3b);
assert(isalmostequal(s3.avg.coh, s3c.avg.coh, 'reltol', 1e-8));
end

% ALSO: check whether it also works with powandcsd data
filename = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg'), 'freq_mtmfft_powandcsd_trl_ctf275');
load(filename);

if 1,
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
s1            = ft_sourceanalysis(cfg, freq);
end

% the following runs through irrespective of the bug being fixed, yet the
% behavior might be unexpected: no single trial output is generated.
if 1,
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.keeptrials = 'yes';
s2             = ft_sourceanalysis(cfg, freq);
end

if 1,
% with the bug unfixed, this gives a problem with the refchan handling, due
% to ft_selectdata only keeping the MEG channels from an early stage
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
s3             = ft_sourceanalysis(cfg, freq);

cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
cfg.keeptrials = 'yes';
s4             = ft_sourceanalysis(cfg, freq);

assert(isalmostequal(s3.avg.coh, s4.avg.coh, 'reltol', 1e-8));
end

% Tthe following cfg should give an error due to the mismatch in
% channels in data and leadfield
if 1,
try,
	cfg           = [];
	cfg.method    = 'dics';
	cfg.frequency = 10;
	cfg.headmodel = vol;
	cfg.grid      = lf;
	s5            = ft_sourceanalysis(cfg, freq);
catch
	% this should be solved with proper error handling
	if ~isempty(strfind(lasterr,'There''s a mismatch'))
		% this is OK
	else
		error('there''s a problem');
	end
end
end

% check whether it also works with pcc as a method
if 1,
cfg           = [];
cfg.method    = 'dics';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
cfg.dics.keepcsd = 'yes';
cfg.dics.realfilter = 'yes';
s3             = ft_sourceanalysis(cfg, freq);

cfg           = [];
cfg.method    = 'pcc';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.grid      = lf;
cfg.channel   = 'MEG'; % this is needed, otherwise there's a detected mismatch in channels vs leadfields
cfg.refchan   = 'BR1';
s3b           = ft_sourceanalysis(cfg, freq);
s3c           = ft_sourcedescriptives([], s3b);
assert(isalmostequal(s3.avg.coh, s3c.avg.coh, 'reltol', 1e-8));
end
