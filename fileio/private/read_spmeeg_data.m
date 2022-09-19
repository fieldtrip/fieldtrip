function dat = read_spmeeg_data(filename, varargin)

% read_spmeeg_data() - import SPM5 and SPM8 meeg datasets
%
% Usage:
%   >> header = read_spmeeg_data(filename, varargin);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'begsample'      first sample to read
%   'endsample'      last sample to read
%   'chanindx'  -    list with channel indices to read
%   'header'    - FILEIO structure header
%
% Outputs:
%   dat    - data over the specified range
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Vladimir Litvak

if nargin < 1
    help read_spmeeg_data;
    return;
end

typenames = {'uint8','int16','int32','float32','float64','int8','uint16','uint32'};
typesizes   = [1  2  4  4 8 1 2 4];

header    = ft_getopt(varargin, 'header');
begsample = ft_getopt(varargin, 'begsample');
endsample = ft_getopt(varargin, 'endsample');
chanindx  = ft_getopt(varargin, 'chanindx');

if isempty(header)
    header = read_spmeeg_header([filename(1:(end-3)) 'mat']);
end

if isempty(begsample), begsample = 1; end
if isempty(endsample), endsample = header.nSamples; end

datatype = 'float32-le';
scale = [];
if isfield(header, 'orig')
    if isfield(header.orig, 'data') && isnumeric(header.orig.data) ...
            && ~isempty(header.orig.data)
        try
            dat = reshape(header.orig.data(chanindx, :, :), length(chanindx), []);
            dat = dat(:, begsample:endsample);
            return;
        end
    end

    if isfield(header.orig, 'datatype')
        datatype = header.orig.datatype;
    elseif isfield(header.orig.data, 'datatype')
        datatype = header.orig.data.datatype;
    end
    if isfield(header.orig, 'scale')
        if isnumeric(header.orig.scale)
            scale = header.orig.scale;
        else
            scale = header.orig.scale.values;
        end
    elseif isfield(header.orig.data, 'scale')
        scale = header.orig.data.scale;
    end
end

stepsize = typesizes(strmatch(strtok(datatype, '-'), typenames));

filename = [filename(1:(end-3)) 'dat'];

fid = fopen_or_error(filename, 'r');
fseek(fid, stepsize*header.nChans*(begsample-1), 'bof');
[dat, siz] = fread(fid, [header.nChans, (endsample-begsample+1)], strtok(datatype, '-'));
fclose(fid);

if ~isempty(chanindx)
    % select the desired channels
    dat = dat(chanindx,:);
end

if ~isempty(scale) && ~ismember(strtok(datatype, '-'), {'float32', 'float64'})
    
    % This is a somewhat complicated mechanism to figure out which scaling
    % coefficients go with which data points in a generic way
    
    trlind = floor(((begsample:endsample)-1)/header.nSamples)+1;

    utrlind = unique(trlind);

    for i = 1:length(utrlind)
        dat(:, trlind == utrlind(i)) = dat(:, trlind == utrlind(i)).* ...
            repmat(squeeze(scale(chanindx,utrlind(i))), 1, sum(trlind == utrlind(i)));
    end

end
