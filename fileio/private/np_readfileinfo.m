function [np_info] = np_readfileinfo (filename, option)
%
% function [np_info] = np_readfileinfo (filename, option)
%
% Purpose:
%
% np_readfileinfo reads out header information from a NEURO PRAX
% data file (*.EEG).
%
% Syntax:
% 
%   [np_info] = np_readfileinfo(filename);
%   [np_info] = np_readfileinfo(filename,'NO_MINMAX');
%
% Input data:
%
%   filename        -   the complete filename with path
%                       (e. g.  C:\Document...\20030716103637.EEG)
%   option          -   if option = 'NO_MINMAX' then physical minima
%   (optional)          and maxima off all channels will not be calculated
%                       (faster for long recordings)
%
% Output data:
%
%   np_info         -   structure
%
%      np_info.filename       -   filename of *.EEG file
%      np_info.pathname       -   pathname of *.EEG file
%      np_info.name           -   patient's name
%      np_info.firstname      -   patient's firstname
%      np_info.birthday       -   patient's birthday
%      np_info.ID             -   identification number
%      np_info.date           -   start date of the recording (format: 'dd-mmm-yyyy')
%      np_info.time           -   start time of the recording (format: 'hh:mm:ss')
%      np_info.duration       -   duration of the recording (global for splitted EEG files!)
%      np_info.setup          -   electrode setup
%      np_info.pmtype         -   type of primary montage
%      np_info.algorithm      -   used algorithm for measurement
%      np_info.channels       -   cell array with channel names
%      np_info.channeltypes   -   cell array with channel types
%      np_info.units          -   cell array with channel units
%      np_info.PhysMin        -   physical minimum for each channel (global for splitted EEG files!)
%      np_info.PhysMax        -   physical maximum for each channel (global for splitted EEG files!)
%      np_info.N              -   number of samples per channel (global for splitted EEG files!)
%      np_info.K              -   number of channels
%      np_info.fa             -   sampling frequency
%      np_info.fp_data        -   filepointer to data
% --- additional SPLITTING information ---
%      np_info.SPLIT_Z        -   number of splitted EEG files (=Z; Z = 1 : no splitting)
%      np_info.SPLIT_filename -   all filenames of splitted EEG file (Zx1 cell-array)
%      np_info.SPLIT_N        -   samples of splitted EEG file (Zx1 array)
%      np_info.SPLIT_fp_data  -   file pointers of splitted EEG file (Zx1 array)
%      np_info.SPLIT_PhysMin  -   physical minimum for each channel and each file (ZxK array)
%      np_info.SPLIT_PhysMax  -   physical maximum for each channel and each file (ZxK array)
%
% Version:  1.3. (2005-09-19)
%
% (1) The field 'units' will be read correctly from the channel header.
%
% (2) The sampling frequency (fa) is identical for all channels and is set
%     to the value of the first channel.
%
% (3) The channel types will be read from the channel header directly.
%
% (4) The field 'time' is set to the correct recording time.
%
% (5) Additional structure fields: 
%     setup      -   the electrode setup (previously "primmon")
%     pmtype     -   the type of the primary montage
%     algorithm  -   the feedback algorithm
%
% (6) No longer available: the structure field "primmon".
%
% (7) Splitted EEG files will be supported.
%
% See also: np_readdata, np_readmarker
%
% eldith GmbH
% Gustav-Kirchhoff-Str. 5
% D-98693 Ilmenau
% Germany
% 02.02.2005


if (1==nargin),
    option='';
end

% -------------------------------------------------------------------------
% Initialisierung der Struktur
% -------------------------------------------------------------------------
np_info=struct('filename','','pathname','','name','','firstname','','birthday','','ID',0,'date','00.00.00',...
               'time','00:00:00','duration','','setup','','pmtype','','algorithm','','channels','','channeltypes','',...
               'units','','PhysMin',0,'PhysMax',0,'N',0,'K',0,'fa',0,'fp_data',0,...
               'SPLIT_Z',1,'SPLIT_filename','','SPLIT_N',[],'SPLIT_fp_data',[],...
               'SPLIT_PhysMin',[],'SPLIT_PhysMax',[]);

% -------------------------------------------------------------------------
% alle zusammengehörenden EEG-Dateien ermitteln (z.B. bei SPLITTING)
% -------------------------------------------------------------------------
[pa,fn,ex]=fileparts(filename);
if strcmp(pa,''),
    pa=pwd;
end
MaskEEGFileNames=[pa filesep fn(1:14) '*.EEG'];
d=dir(MaskEEGFileNames);
e=cell(length(d),1);
for i=1:length(d)
    e{i}=[pa filesep d(i).name];
end
filenames=sort(e);

% -------------------------------------------------------------------------
% SPLIT-Procedures
% -------------------------------------------------------------------------
np_info.SPLIT_Z=length(filenames);
np_info.SPLIT_filename=cell(np_info.SPLIT_Z,1);
np_info.SPLIT_N=zeros(np_info.SPLIT_Z,1);
np_info.SPLIT_fp_data=zeros(np_info.SPLIT_Z,1);

for Z=1:np_info.SPLIT_Z
    % Infos für nächste EEG-Datei lesen
    np_splitted_info=np_read_splitted_fileinfo(filenames{Z},option);
    
    % fixe Daten in np_info eintragen
    if (1==Z),
        np_info.filename=np_splitted_info.filename;
        np_info.pathname=np_splitted_info.pathname;
        np_info.name=np_splitted_info.name;
        np_info.firstname=np_splitted_info.firstname;
        np_info.birthday=np_splitted_info.birthday;
        np_info.ID=np_splitted_info.ID;
        np_info.date=np_splitted_info.date;
        np_info.time=np_splitted_info.time;
        np_info.setup=np_splitted_info.setup;
        np_info.pmtype=np_splitted_info.pmtype;
        np_info.algorithm=np_splitted_info.algorithm;
        np_info.K=np_splitted_info.K;
        np_info.fa=np_splitted_info.fa;
        np_info.fp_data=np_splitted_info.fp_data;
        np_info.channels=np_splitted_info.channels;
        np_info.channeltypes=np_splitted_info.channeltypes;
        np_info.units=np_splitted_info.units;
        np_info.PhysMin=zeros(1,np_splitted_info.K);
        np_info.PhysMax=zeros(1,np_splitted_info.K);
        np_info.SPLIT_PhysMin=zeros(np_info.SPLIT_Z,np_splitted_info.K);
        np_info.SPLIT_PhysMax=zeros(np_info.SPLIT_Z,np_splitted_info.K);
    end
    
    % SPLIT-Infos zusammenstellen
    np_info.SPLIT_filename{Z}=np_splitted_info.filename;
    np_info.SPLIT_N(Z)=np_splitted_info.N;
    np_info.SPLIT_fp_data(Z)=np_splitted_info.fp_data;
    np_info.SPLIT_PhysMin(Z,:)=np_splitted_info.PhysMin;
    np_info.SPLIT_PhysMax(Z,:)=np_splitted_info.PhysMax;

    % Samplezähler erhöhen
    np_info.N=np_info.N+np_info.SPLIT_N(Z);
    
end

% global duration berechnen
np_info.fa;
np_info.duration=num2str(np_info.N./np_info.fa);

% physikalische Minima und Maxima berechnen
if (1~=strcmp(upper(option),'NO_MINMAX')),
    for k=1:np_info.K
        np_info.PhysMin(k)=+1e38;
        np_info.PhysMax(k)=-1e38;
        for Z=1:np_info.SPLIT_Z
            np_info.PhysMin(k)=min([np_info.PhysMin(k) np_info.SPLIT_PhysMin(Z,k)]);
            np_info.PhysMax(k)=max([np_info.PhysMax(k) np_info.SPLIT_PhysMax(Z,k)]);
        end
    end
end
