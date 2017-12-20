function [dataout] = ft_headposavg(cfg, datain)

% Use as
%   [trans] = ft_headposavg(cfg)
% Read continous head position estimation from fif file and return average
% 4x4 transformation matrix from dewar to head coordinates, based on the the
% continous head position. cfg is a configuration structure that should 
% contain the filename of a fif file processed with Neuromag MaxFilter.
%
% The configuration should contain:
%  cfg.dataset        = a string (e.g. 'mydata.fif')
%
% The configuration can optionally contain:
%  cfg.writetransfif  = Write a fif file containing the averaged
%                       transformation from dewar to head coordinates,
%                       'yes'/'no' (default = 'no')
%  cfg.outfname       = string, filename of output head transformation fif
%                       file (default = '[dataset]-trans.fif'
%  cfg.outdir         = string with output directory (NB! default = current
%                       dir)
%  cfg.summarystat    = 'mean' or 'median', how to summarize headpos
%                       (default = 'mean')
%  cfg.visualize      = string, 'yes'/'no' plot traces of the head position and 
%                       rotation (default = 'no')
%  cfg.saveplot       = string, 'yes' or 'no' to save the visualization. Will
%                       save to same folder as cfg.outdir (default = 'yes')
%  
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_GETCONTHEADTRANS

% Copyrights (C) 2017, Mikkel C. Vinding (mikkel.vinding@ki.se).
%
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
datain = ft_checkdata(datain, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'deprecated',  {'normalizecov', 'normalizevar'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% ensure that the required options are present and valid
cfg = ft_checkconfig(cfg, 'required', {'dataset'});
cfg = ft_checkopt(cfg, 'dataset', 'char');

infpath    = ft_getopt(cfg, 'dataset');       
F = strsplit(infpath, '/') ;
inFname = F{end};

try
    cfg = ft_checkconfig(cfg, 'required', {'outfname'});
    outFname = ft_getopt(cfg, 'outfname');
catch
    O = strsplit(inFname,'.');
    outFname = O{1};
    outFname = [outFname,'-trans.fif'];
end

% ensure that the options are valid
% cfg = ft_checkopt(cfg, 'writetransfif', 'char', {'yes', 'no'});
% cfg = ft_checkopt(cfg, 'dataset', 'char');

% set the defaults 
cfg.writetransfif   = ft_getopt(cfg, 'writetransfif',   'no');
cfg.outfname        = ft_getopt(cfg, 'outfname',   outFname);
cfg.outdir          = ft_getopt(cfg, 'outdir',   pwd);
cfg.summarystat     = ft_getopt(cfg, 'summarystat',   'mean');
cfg.visualize       = ft_getopt(cfg, 'visualize',   'no');
cfg.saveplot        = ft_getopt(cfg, 'saveplot',   'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'writetransfif', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'summarystat', 'char', {'mean', 'median'});

% get the options
method    = ft_getopt(cfg, 'summarystat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfgg = [];
cfgg.dataset    = cfg.dataset;
cfgg.returndata = 'mat';
[H, t, fs] = ft_getcontheadtrans(cfgg);

% Convert from mat to angles
A = zeros(length(H),3);
for i = 1:length(H)
    A(i,:) = rot2ang(H(:,:,i));
end

% Calculate summary
if strcmp(method,'mean')
    A_mean = mean(A);
    rot3d = ang2rot(A_mean(1),A_mean(2),A_mean(3));
    T = mean(H,3);
    T(1:3,1:3) = rot3d; 
elseif strcmp(method, 'median')
    A_medi = median(A);
    rot3d = ang2rot(A_medi(1),A_medi(2),A_medi(3));
    T = median(H,3);
    T(1:3,1:3) = rot3d;
else
    error('Something is wrong!');
end

dataout = [T];

% Write fif file
if strcmp(writetransfif, 'yes')
    fiffname = fullfile(cfg.outdir,cfg.outfname);
    fid = fiff_start_file(fiffname);
    fiff_start_file(fid,T)
    fiff_end_file(fid)
    fprintf('Wrote transfomation to %s\n',fiffname)
end

% Plots
if strcmp(cfg.visualize, 'yes')
    
    % Read headpos
    head_dev = ft_read_headshape(cfg.dataset,'coordsys','dewar');
    fidOrigDev = head_dev.fid.pos;

    fid1 = zeros(3,length(H));
    fid2 = zeros(3,length(H));
    fid3 = zeros(3,length(H));

% Must have a way to indetify NAS LPA RPA here...       [!] why?

    % [Horig] = ft_headcoordinates(fidOrigHead(2,:),fidOrigHead(1,:),fidOrigHead(3,:),head_dev.coordsys);

    tic
    for i = 1:length(H)
        fid1(:,i) = ft_warp_apply(H(:,:,i),fidOrigDev(1,:)); 
        fid2(:,i) = ft_warp_apply(H(:,:,i),fidOrigDev(2,:)); 
        fid3(:,i) = ft_warp_apply(H(:,:,i),fidOrigDev(3,:)); 
    end
    
    fid1avg = ft_warp_apply(T,fidOrigDev(1,:));
    fid2avg = ft_warp_apply(T,fidOrigDev(2,:));
    fid3avg = ft_warp_apply(T,fidOrigDev(3,:));
    
    cc = circumcenter(fid1, fid2, fid3);
    cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';
    toc

    ccAvg = circumcenter(fid1avg', fid2avg', fid3avg');
    ccAvg_rel = [ccAvg - cc(:,1)]';

    
    % plot translations
    figure();
    subplot(2,1,1); plot(t, cc_rel(:,1:3)*1000) % in mm
    title('Position'); xlabel('Time (s)'); ylabel('mm');
    hold on
    h1 = refline([0 ccAvg_rel(1)*1000]);
    h1.Color = 'b';
    h2 = refline([0 ccAvg_rel(2)*1000]); 
    h2.Color = 'r';
    h3 = refline([0 ccAvg_rel(3)*1000]);
    h3.Color = 'y';

    % plot rotations
    subplot(2,1,2); plot(t, cc_rel(:,4:6))
    title('Roation'); xlabel('Time (s)'); ylabel('rad'); % [!] I think this is radians, but I am not sure

    % Save plot
    if strcmp(cfg.saveplot,'yes')
        kk = strsplit(inFname,'.');
        plotfname = kk{1};
        set(gcf, 'PaperType', 'a4');
        
        print(gcf, '-dpng', strcat(cfg.outdir,'/',plotfname,'.png'));
        orient landscape;
        print(gcf, '-dpdf', strcat(cfg.outdir,'/',plotfname,'.pdf'));
      close
    end
    
end

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
