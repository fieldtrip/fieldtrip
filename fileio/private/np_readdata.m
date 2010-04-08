function [np_data] = np_readdata (filename, np_info, idx_begin, data_length, option)
%
% function [np_data] = np_readdata (filename, idx_begin, data_length, option)
%
% np_readdata reads data from a NEURO PRAX data file (*.EEG).
%
%
% Syntax:
%
%   [np_data] = np_readdata(filename,idx_begin,data_length,'samples');
%   [np_data] = np_readdata(filename,idx_begin,data_length,'time');
%
% Input data:
%
%   filename    -   the complete filename with path
%                   (e. g.  C:\Document...\20030716103637.EEG)
%   idx_begin   -   the start index of the data block to be read
%   data_length -   the length of the data block to be read
%   option      -   if option = 'samples':
%                       the data block starts at sample 'idx_begin' from the recording;
%                       data_length is the number of samples to be read
%                   if option = 'time':
%                       the data block starts at time 'idx_begin' from the recording;
%                       data_length is the number of seconds to be read
%
%                   To read all data use: idx_start = 0, data_length = inf, option =
%                   'samples'.
%
% Output data:
%
%   np_data     -   structure
%      np_data.data     -   data matrix of unipolar raw data
%                           dimension of the matrix: (NxK)
%                           N: number of samples
%                           K: number of channels (each column is one channel)
%      np_data.t        -   discrete time vector for the recording
%
% Version:  1.2. (2005-02-02)
%
% See also: np_readfileinfo, np_readmarker
%
% eldith GmbH
% Gustav-Kirchhoff-Str. 5
% D-98693 Ilmenau
% Germany
% 02.02.2005

% Daten initialisieren, Startsamples und Blockl‰nge ermitteln
if (0==idx_begin) && (inf==data_length)
    N_start=0;
    N = np_info.N;                 
else
    switch upper(option)
        case 'SAMPLES', 
            N_start = idx_begin;
            N       = data_length;
        case 'TIME',
            N_start = round(idx_begin*np_info.fa);
            N       = round(data_length*np_info.fa);
        otherwise,
            error ('Bad specification for ''option''');
    end
    if (N_start<0)
        error ('idx_begin is to small.');
    end
    if ((N_start+N-1)>(np_info.N-1))
        error('data_length is to big.');
    end
end
np_data.data=zeros(N,np_info.K);
np_data.t=(N_start:N_start+N-1)'./np_info.fa;

% -------------------------------------------------------------------------
% Messadten einlesen
%
% sequentielles Datenformat in Datei (K-Kanalindex, N-Sampleindex):
% x11 x21 x31 ... xK1 ...  
% x12 x22 x32 ... xK2 ...
% ...
% x1(N-1) x2(N-1) ... xK(N-1) ...
% x1N x2N x3N ... xKN
%
% einlesen in eine Matrix mit der Dimension: (KxN)
% anschlieﬂend transponieren in eine Matrix der Dimension: (NxK)
%
% -------------------------------------------------------------------------
fid=fopen(filename,'r');
if fid==-1,
    error(['Unable to open file "' filename '" . Error code: ' ferror(fid)]);
end
status=fseek(fid,np_info.fp_data+N_start*np_info.K*4,'bof');    % 4 Byte pro Sample
if status~=0,
    fclose(fid);
    error('Unable to set filepointer to begin of data block.');
end
[np_data.data,count]=fread(fid,[np_info.K N],'float');      % lese Messdaten


%if count~=N*np_info.K,
%    fclose(fid);
%    error('Number of read samples unequal to product N*K.');
%end
np_data.data=np_data.data';         % transponieren
fclose(fid);
