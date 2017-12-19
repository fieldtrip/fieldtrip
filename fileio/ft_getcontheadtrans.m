function [varargout] = ft_getcontheadtrans(cfg)

% FT_GETCONTHEADTRANS imports transformation of head postion read from the
% countinous head position estimation from Neuromag data on which MaxFilter
% has been performed. This can be from a file where headpostion, but not 
% head movement compensation, has been performed ('-quat' file) or any file
% on which SSS/tSSS with head movement compensation has been performed.
%
% Use as
%   [transform, time, sfreq] = ft_getcontheadtrans(cfg)
% where cfg is a configuration structure that must contain the filename of 
% a file on which Neuromag MaxFilter has been performed. Will return a
% 4x4xN transformation matrix or 6xN quaternions describing the transform 
% from dewar coordinates to head coordinates, where N is the total number
% of samples.
%
%  cfg.dataset     = filename as a string (e.g. 'data-tsss.fif')
%  cfg.returndata  = 'quat' or 'mat', should the function return quaternion
%                   (6xN matrix) or transformation (4x4xN) matrix 
%                   (default = 'mat')
%
% The configuration can optionally contain
%  cfg.resample   = 'yes'/'no', downsample data (default = 'no')
%  cfg.resamplefs = frequency at which the data will be resampled (default = 100 Hz)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING, FT_READ_HEADSHAPE

% Copyright (C) 2017, Mikkel C. Vinding (mikkel.vinding@ki.se)
%
% Here comes the Revision tag, which is auto-updated by the version control system
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar    datain % this reads the input data in case the user specified the cfg.inputfile option
ft_preamble provenance datain % this records the time and memory usage at the beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
% datain = ft_checkdata(datain, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'deprecated',  {'normalizecov', 'normalizevar'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% ensure that the required options are present
%cfg = ft_checkconfig(cfg, 'required', {'method', 'foi', 'tapsmofrq'});

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'dataset', 'char');
cfg = ft_checkopt(cfg, 'returndata', 'char', {'mat', 'quat'});

% get the options
fname    = ft_getopt(cfg, 'dataset');        % there is no default

% set the defaults
cfg.resample    = ft_getopt(cfg, 'resample',   'no');
cfg.resamplefs  = ft_getopt(cfg, 'resamplefs',   100);
cfg.returndata  = ft_getopt(cfg, 'returndata',   'mat');

%% Read data
hdr = ft_read_header(fname);
if any(strfind(fname, 'quat'))
    quatChans = find(~cellfun(@isempty, strfind(hdr.label, 'QUAT')));
    gofChans = find(~cellfun(@isempty, strfind(hdr.label, 'QUAT007')));    
else %if any(strfind(fname, 'tsss'))
    quatChans = find(~cellfun(@isempty, strfind(hdr.label, 'CHPI')));
    gofChans = find(~cellfun(@isempty, strfind(hdr.label, 'CHPI007')));
    if ~isempty(quatChans)
        fprintf('Found CHPIXXXX channels in %s\n', fname);
    else
        error('Error: Cannot find CHPIXXXX channels in %s\n', fname);
    end
end

gof = ft_read_data(fname, 'chanindx', gofChans);

cfgd = [];
cfgd.dataset    = fname;
cfgd.channel     = quatChans(1:6); %Only first 6 channels are quaternions
cfgd.continous   = 'yes';
dat = ft_preprocessing(cfgd);

% In case cHPI is not started first 
if any(gof < 0.98)                  % cutoff for MaxFilter to work is 0.98
    begsample = find(gof > 0.98, 1);          
    cfgs = [];
    cfgs.begsample = begsample;
    cfgs.endsample = dat.sampleinfo(2);
    
    dat = ft_redefinetrial(cfgs,dat);
end

if strcmpi(cfg.resample, 'yes')
    cfgrs = [];
    cfgrs.resamplefs = cfg.resamplefs;
    dat = ft_resampledata(cfg,dat);
end

q = dat.trial{:};
t = dat.time{:};
fs = dat.fsample;

if strcmp(cfg.returndata, 'quat')
    outdat = q;
elseif strcmp(cfg.returndata, 'mat')        % Make transformation

    H = zeros(4,4,length(q));
    for i = 1:length(q)
        H(:,:,i) = quaternion(q(:,i));
    end
    
    outdat = H;
else
    error('Something has gone wrong!')
end

% 
% % Read headpos
% head = ft_read_headshape(fname,'coordsys','head');
% fidOrigHead = head.fid.pos;
% head_dev = ft_read_headshape(fname,'coordsys','dewar');
% fidOrigDev = head_dev.fid.pos;
% 
% fid1 = zeros(3,length(H));
% fid2 = zeros(3,length(H));
% fid3 = zeros(3,length(H));
% 
% % Must have a way to indetify NAS LPA RPA here...       [!]
% 
% [Horig] = ft_headcoordinates(fidOrigHead(2,:),fidOrigHead(1,:),fidOrigHead(3,:),head_dev.coordsys);
% 
% tic
% for i = 1:length(H)
%     fid1(:,i) = ft_warp_apply(H(:,:,i),fidOrigDev(1,:)); 
%     fid2(:,i) = ft_warp_apply(H(:,:,i),fidOrigDev(2,:)); 
%     fid3(:,i) = ft_warp_apply(H(:,:,i),fidOrigDev(3,:)); 
% end
% toc
% 
% figure; ft_plot_headshape(head_dev); hold on
% plot3(fid1(1,:),fid1(2,:),fid1(3,:),'xk'); hold on
% plot3(fid2(1,:),fid2(2,:),fid2(3,:),'xk'); 
% plot3(fid3(1,:),fid1(2,:),fid3(3,:),'xk'); hold off
% 
% cc = circumcenter(fid1, fid2, fid3);
% 
% cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';
%  
% % plot translations
% figure(); 
% subplot(2,1,1); plot(t, cc_rel(:,1:3)*1000) % in mm
% title('Position'); xlabel('Time (s)'); ylabel('mm');
% 
% % plot rotations
% subplot(2,1,2); plot(t, cc_rel(:,4:6))
% title('Roation'); xlabel('Time (s)'); ylabel('rad'); % [!] I think this is radians, but I am not sure
% 
% 
% % do your stuff...
% dataout = [];

% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function

% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig         % this converts the config object back into a struct and can report on the unused fields
ft_postamble previous   datain   % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble provenance dataout  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble history    dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar    dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option

%% VARARGOUT
if nargout>0
  mOutputArgs{1} = outdat;
  mOutputArgs{2} = t;
  mOutputArgs{3} = fs;
  [varargout{1:nargout}] = mOutputArgs{:};
  clearvars -except varargout
else
  clear
end

