function [montage, cfg] = ft_prepare_ODtransformation(cfg, data)

% FT_PREPARE_ODTRANSFORMATION returns the transformation matrix from 
% optical densities (OD) to chromophore concentrations such as (de-)
% oxygenated hemoglobin.
%
% Use as
%   [montage] = ft_prepare_ODtransformation(cfg, data);
%
% It is neccessary to input the data on which you want to perform the
% inverse computations, since that data generally contain the optode
% information and information about the channels that should be included in
% the transformation. The data structure can be either obtained
% from FT_PREPROCESSING, FT_FREQANALYSIS, FT_TIMELOCKANALYSIS or
% FT_COMPONENTANALYSIS. If the data is empty, all channels will be included
% in transformation.
%
% The configuration should contain
%  cfg.channel            = Nx1 cell-array with selection of channels
%                           (default = 'nirs'), see FT_CHANNELSELECTION for
%                           more details
%
% Optional configuration settings are
%  cfg.age                = scalar, age of the subject (necessary to
%                           automatically select the appropriate DPF, or
%  cfg.dpf                = scalar, differential path length factor
%  cfg.dpffile            = string, location to a lookup table for the
%                           relation between participant age and DPF
%                           
%
%
% Note that the DPF might be different across channels, and is usually
% contained in the optode structure contained in the data.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_TRANSFORM_ODS, FT_APPLY_MONTAGE

% options to be implemented:
% 
% The NIRS positions can be present in the data or can be specified as
%   cfg.opto          = structure with optode positions, see FT_DATATYPE_SENS
%   cfg.optofile      = name of file containing the optode positions, see FT_READ_SENS
%
% cfg.siunits  = ft_getopt(cfg, 'siunits', 'no');   % yes/no, ensure that SI units are used consistently
%
% cfg.logarithm           = string, can be 'natural' or 'base10' (default =

% You are using the FieldTrip NIRS toolbox developed and maintained by 
% Artinis Medical Systems (http://www.artinis.com). For more information
% on FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% 
% This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% 
% Creative Commons Attribution-ShareAlike 4.0 International License:
% -----------------------------------
% You are free to:
% 
%     Share � copy and redistribute the material in any medium or format
%     Adapt � remix, transform, and build upon the material
%     for any purpose, even commercially.
% 
%     The licensor cannot revoke these freedoms as long as you follow the 
%     license terms.
% 
% Under the following terms:
% 
%     Attribution � You must give appropriate credit, provide a link to 
%                    the license, and indicate if changes were made. You 
%                    may do so in any reasonable manner, but not in any way 
%                    that suggests the licensor endorses you or your use.
% 
%     ShareAlike � If you remix, transform, or build upon the material, 
%                   you must distribute your contributions under the same 
%                   license as the original.
% 
%     No additional restrictions � You may not apply legal terms or 
%                                   technological measures that legally 
%                                   restrict others from doing anything the 
%                                   license permits.
% 
% -----------------------------------
% 
% This toolbox is not to be used for medical or clinical purposes.
% 
% Copyright (c) 2015 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmer: J�rn M. Horschig
% $Id: ft_prepare_ODtransformation.m 10169 2015-02-05 12:56:47Z jorhor $

revision = '$Id: ft_prepare_ODtransformation.m 10169 2015-02-05 12:56:47Z jorhor $';


% check if the input data is valid for this function
data = ft_checkdata(data);

% set the defaults
cfg.channel = ft_getopt(cfg, 'channel', 'nirs');
cfg.age     = ft_getopt(cfg, 'age', []);
cfg.dpf     = ft_getopt(cfg, 'dpf', []);

% get the optode definition
%sens = ft_fetch_sens(cfg, data);
if ~isfield(data, 'opto')
  error('no optode structure found in the data');
end

sens = data.opto;

% select the appropriate channels
if isfield(data, 'topolabel')
  % the data reflects a componentanalysis, where the topographic and the
  % timecourse labels are different
  cfg.channel = ft_channelselection(cfg.channel, data.topolabel);
elseif isfield(data, 'label')
  % In the subsequent code, the matching channels in the sensor array and
  % in the configuration will be selected. To ensure that these channels
  % are also present in the data, update the configuration to match the data.
  cfg.channel = ft_channelselection(cfg.channel, data.label);
else
  % update the selected channels based on the optode definition
  cfg.channel = ft_channelselection(cfg.channel, sens.label);
end
cfg.channel = ft_channelselection('nirs', cfg.channel); % only NIRS valid

% check current data type by looking what's inside of []
typeidx = regexp(cfg.channel, sprintf('%s%s', regexptranslate('wildcard', '[*]'), '$'));

% remove channels that do not obey to labelling guidelines
illegalidx = cellfun(@isempty, typeidx);
cfg.channel(illegalidx) = [];

if isempty(cfg.channel)
  warning('no valid NIRS channels found')
  return;
end

% channel indices wrt optode structure
chanidx     = match_str(sens.label, cfg.channel);

% check on DPF values
if ~isempty(cfg.age) && ~isempty(cfg.dpf)
  error('unsure if age or dpf in configuration should take precedence');
elseif ~isempty(cfg.age)
  error('the use of cfg.age is not implemented, yet');
elseif ~isempty(cfg.dpf)
  dpfs = repmat(cfg.dpf, size(cfg.channel));
else
  dpfs = sens.DPF(chanidx);
end

% which chromophores are desired
chromophoreIdx = [3 2];
chromophoreName = {'O2Hb' 'HHb'};

%% transform to concentrations or to OD
% read in the coefficient table
% TODO FIXME put this into an own function to read this out
fid = fopen(fullfile(fileparts(mfilename('fullpath')), 'private', 'Cope_ext_coeff_table.txt'));
coefs = cell2mat(textscan(fid, '%f %f %f %f %f'));

% extract all transceivers that are relevant here
transceivers   = sens.transceiver(chanidx, :);
transmitteridx = transceivers>0;
receiveridx    = transceivers<0;
fiberidx       = transmitteridx | receiveridx;

% extract the wavelengths
wavelengths  = sens.wavelength(transceivers(transmitteridx));  
wlidx = bsxfun(@minus, coefs(:, 1), wavelengths);

% find the relevant channel combinations
chancmb = logical([]); % we need to loop here to get size of chancmb
chanUsed = zeros(numel(chanidx), 1);
for c=1:numel(chanidx)
  % skip channel if it already was in another channel combination 
  if chanUsed(c)
    continue;
  end
  % compute the channel combinations
  tupletidx = sum(bsxfun(@minus, fiberidx, fiberidx(c, :))~=0, 2)==0;
  chanUsed = chanUsed|tupletidx;
  chancmb(:, end+1) = tupletidx;
end

% do the transformation
tra      = zeros(size(chancmb, 2)*numel(chromophoreIdx), size(cfg.channel, 1));
labelnew = cell(1, size(chancmb, 2)*numel(chromophoreIdx));

% transformation has to be done per channel combination
chanidx = []; % we will use this variable here
for c=1:size(chancmb, 2)
  
  % find all other channels of the same transmitter / receiver pair
  chanidx = chancmb(:, c);
  
  % extract coefficient idx of these channels
  [coefidx, colidx] = find(wlidx(:, chanidx)==0);
  
  % compute the transmitter/receiver distance in cm
  dist = sqrt(sum(diff(sens.fiberpos(fiberidx(c, :), :)).^2));
  
  % select dpf
  dpf = mean(dpfs(chanidx));
  
  % compute distance factor for correct scaling
  distfactor = dist/10*dpf/10; % in mm
  
  % select the respective coefficients
  eps=coefs(coefidx,chromophoreIdx)/10; % in mm
  
  % add to transformation matrix
  traIdx = 1+(c-1)*numel(chromophoreIdx):(c)*numel(chromophoreIdx);
  tra(traIdx, chanidx) = pinv(eps)/(distfactor);
  
  % create the new labels for the channels
  chanorg = cfg.channel(find(chanidx, 1, 'first'));
  for n=1:numel(chromophoreIdx)
    labelnew(traIdx(n), 1) = regexprep(chanorg, sprintf('%s%s', regexptranslate('wildcard', '[*]'), '$'), sprintf('[%s]', chromophoreName{n}));
  end
  
  
end

%% create output
montage = [];
montage.labelorg = cfg.channel;
montage.labelnew = labelnew;
montage.tra      = tra;

