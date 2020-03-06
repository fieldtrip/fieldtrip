function [dat, ref] = ft_preproc_rereference(dat, refchan, method, handlenan, leadfield)

% FT_PREPROC_REREFERENCE computes the average reference over all EEG channels
% or rereferences the data to the selected channels
%
% Use as
%   [dat] = ft_preproc_rereference(dat, refchan, method, handlenan,leadfield)
% REST example: [dat] = ft_preproc_rereference(dat, refchan,
%                       'rest',[],leadfield);
% where
%   dat        data matrix (Nchans X Ntime)
%   refchan    vector with indices of the new reference channels, or 'all'
%   method     string, can be 'avg', 'median', or 'rest'
%              if select 'rest','leadfield' is required.
%              The leadfield can be a matrix (channels X sources)
%              which is calculated by using the forward theory, based on
%              the electrode montage, head model and equivalent source
%              model. It can also be the output of ft_prepare_leadfield.m
%              (e.g. lf.leadfield) based on real head modal using FieldTrip.
%   handlenan  boolean, can be true or false
%
% If the new reference channel is not specified, the data will be
% rereferenced to the average of all channels.
%
% See also PREPROC

% Copyright (C) 1998-2017, Robert Oostenveld
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

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% determine the new reference channels
if nargin<2 || isempty(refchan) || (ischar(refchan) && strcmp(refchan, 'all'))
    refchan = 1:Nchans;
end

if nargin<3 || isempty(method)
    method = 'avg';
end

if nargin<4 || isempty(handlenan)
    handlenan = false;
end

if nargin < 5 && isequal(method,'rest')
    ft_error('Leadfield is required to re-refer to REST');
elseif nargin == 5 && isempty(leadfield)
    ft_error('Leadfield is empty');
end

hasnan = any(any(isnan(dat(refchan,:))));

if hasnan && handlenan
    % preprocessing works differently if channels contain NaN
    switch method
        case 'avg'
            ref = nanmean(dat(refchan,:), 1);
        case 'median'
            ref = nanmedian(dat(refchan,:), 1);
        case 'rest'
            ft_error('channels contain NaN, and REST method is not supported');
        otherwise
            ft_error('unsupported method')
    end % switch
else
    % preprocessing fails on channels that contain NaN
    if any(isnan(dat(:)))
        ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
    end
    % compute the average value over the reference channels
    % or re-referencing to REST
    switch method
        case 'avg'
            ref = mean(dat(refchan,:), 1);
        case 'median'
            ref = median(dat(refchan,:), 1);
        case 'rest' % re-referencing to REST
            % get the leafield matrix
            if isnumeric(leadfield)
                G = leadfield';
            elseif isstruct(leadfield)
                try
                    Npos = size(leadfield.pos,1);
                    m = 1;
                    for i = 1:Npos
                        if ~isempty(leadfield.leadfield{1,i})
                            lf_X(:,m) = leadfield.leadfield{1,i}(:,1); % X orientation of the dipole.
                            lf_Y(:,m) = leadfield.leadfield{1,i}(:,2); % Y orientation of the dipole.
                            lf_Z(:,m) = leadfield.leadfield{1,i}(:,3); % Z orientation of the dipole.
                            m = m + 1;
                        end
                    end
                    G = [lf_X,lf_Y,lf_Z]';
                    % the leadfield matrix (sources*3 X chans), which
                    % contains the potential or field distributions on all
                    % sensors for the x,y,z-orientations of the dipole.
                catch
                    try
                        Npos = size(leadfield.pos,1);
                        m = 1;
                        for i = 1:Npos
                            if ~isempty(leadfield.leadfield{1,i})
                                G(:,m) = leadfield.leadfield{1,i};% potential of the dipole.
                                m = m + 1;
                            end
                        end
                        G = G';
                    catch
                        ft_error('leadfiled is not calculated by ''ft_prepare_leadfield.m'' (dipoles contain x,y,z-orientations)?');
                    end
                    
                end
            elseif iscell(leadfield)
                try
                    Npos = length(leadfield);
                    m = 1;
                    for i = 1:Npos
                        if ~isempty(leadfield{1,i})
                            lf_X(:,m) = leadfield{1,i}(:,1); % X orientation of the dipole.
                            lf_Y(:,m) = leadfield{1,i}(:,2); % Y orientation of the dipole.
                            lf_Z(:,m) = leadfield{1,i}(:,3); % Z orientation of the dipole.
                            m = m + 1;
                        end
                    end
                    G = [lf_X,lf_Y,lf_Z]';
                    % the leadfield matrix (sources*3 X chans), which
                    % contains the potential or field distributions on all
                    % sensors for the x,y,z-orientations of the dipole.
                catch
                    try
                        Npos = size(leadfield.pos,1);
                        m = 1;
                        for i = 1:Npos
                            if ~isempty(leadfield{1,i})
                                G(:,m) = leadfield{1,i};% potential of the dipole.
                                m = m + 1;
                            end
                        end
                        G = G';
                    catch
                        ft_error('leadfiled is not calculated by ''ft_prepare_leadfield.m'' (dipoles contain x,y,z-orientations)?');
                    end
                end
            end
            temp_dat = dat(refchan,:);
            dat = rest_refer(temp_dat,G(:,refchan)); % re-rerefencing to REST
            ref = [];
        otherwise
            ft_error('unsupported method')
    end % switch
end

if ~isequal(method,'rest')
    % subtract the new reference from the data
    for chan=1:Nchans
        dat(chan,:) = dat(chan,:) - ref;
    end
end
% -------------------------------------------------------------------------
% subfunction
function [data_z] = rest_refer(data,G)
%   Main function of Reference Electrode Standardization Technique
%   Input:
%         data:  The EEG potentials,channels X time points,
%                e.g. 62 channels X 10000 time points.
%         G: Lead Field matrix, sources X channels, e.g. 3000 sources X 62 channels.
%   Output:
%         data_z: The EEG potentials with zero reference,
%                channels X time points.
%
%  Edit by Li Dong (Oct. 26, 2017)
%   For more see http://www.neuro.uestc.edu.cn/rest/
%   Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
%                       Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305

if nargin < 2
    error('Please input the Lead Field matrix!');
end
Nchans = size(data,1);
data = data - repmat(mean(data),Nchans,1); % average

G = G';
if size(data,1)~=size(G,1)
    error('No. of channels of leadfield matrix and data are NOT equal!');
end

Gar = G - repmat(mean(G),size(G,1),1);
data_z = G * pinv(Gar,0.05) * data;  % the value 0.05 is for real data;for
% simulated data, it may be set as zero.

data_z = data + repmat(mean(data_z),size(G,1),1); % V = V_avg + AVG(V_0)