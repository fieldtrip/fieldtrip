function [dat] = read_nmc_archive_k_data(datafile, hdr, begsample, endsample, channelsel)

% READ_NMC_ARCHIVE_K_DATA reads data from nmc_archive_k datasets
%
% Used in read_data as
%   dat = read_nmc_archive_k_data(datafile, hdr, begsample, endsample, channelsel);
%
%
% This function specifically only reads data from one of the archived
% datasets of the Neurophysiological Mechanisms of Cognition group of
% Eric Maris, at the Donders Centre for Cognition, Radboud University,
% Nijmegen, the Netherlands. It should not be used for any other data
% format.
%
%
%

% Copyright (C) 2009, Roemer van der Meij
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



% Sanity checks (to be added)



% Getting data-directory out of data-filename and check
datadir = datafile(1:(strfind(datafile, 'eeg.noreref')+11));
if exist(datadir,'dir') ~= 7
    ft_error('no proper data-directory provided, please check your paths');
end

% Getting session name + path out of data-filename
sessionpath = datafile(1:(end-4));


% Build array containing channel file extensions to be used in reading the data
channelext = [];
% Remove the 'CH' from channel labels
chanlabnum = hdr.label;
for ichan = 1:length(chanlabnum)
    chanlabnum{ichan} = chanlabnum(3:end);
end
for ichan = 1:length(channelsel)
    if length(hdr.label{channelsel(ichan)}) == 3
        channelext{ichan} = ['.00' hdr.label{channelsel(ichan)}(3:end)];
    elseif length(hdr.label{channelsel(ichan)}) == 4
        channelext{ichan} = ['.0' hdr.label{channelsel(ichan)}(3:end)];
    elseif  length(hdr.label{channelsel(ichan)}) == 5
        channelext{ichan} = ['.' hdr.label{channelsel(ichan)}(3:end)];
    end
end

% Loop over channels while reading in data from file, 1 file per channel
dat = zeros(length(channelsel),(endsample-begsample+1));
for ichan = 1:length(channelsel)
    channelfile = [sessionpath channelext{ichan}];
    datafid = fopen_or_error(channelfile,'r','l');
    fseek(datafid,(hdr.nBytes*(begsample-1)),'bof');
    dat(ichan,:) = fread(datafid,(endsample-begsample+1),hdr.dataformat)';
    fclose(datafid);
end
