% openbdf() - Opens an BDF File (European Data Format for Biosignals) in MATLAB (R)
%
% Usage:
%   >> EDF=openedf(FILENAME)
%
% Note: About EDF -> www.biosemi.com/faq/file_format.htm
%
% Author: Alois Schloegl, 5.Nov.1998
%
% See also: readedf()

%   Copyright (C) 1997-1998 by Alois Schloegl
%   a.schloegl@ieee.org
%   Ver 2.20    18.Aug.1998
%   Ver 2.21    10.Oct.1998
%   Ver 2.30     5.Nov.1998
%
%   For use under Octave define the following function
% function s=upper(s); s=toupper(s); end;

% V2.12    Warning for missing Header information
% V2.20    EDF.AS.* changed
% V2.30    EDF.T0 made Y2K compatible until Year 2090

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Name changed for bdf files Sept 6,2002  T.S. Lorig
% Header updated for EEGLAB format (update web link too) - Arnaud Delorme 14 Oct 2002

function [DAT,H1]=openbdf(FILENAME)

SLASH='/';      % defines Seperator for Subdirectories
BSLASH=char(92);

cname=computer;
if cname(1:2)=='PC' SLASH=BSLASH; end;

fid=fopen_or_error(FILENAME,'r','ieee-le');          

EDF.FILE.FID=fid;
EDF.FILE.OPEN = 1;
EDF.FileName = FILENAME;

PPos=min([max(find(FILENAME=='.')) length(FILENAME)+1]);
SPos=max([0 find((FILENAME=='/') | (FILENAME==BSLASH))]);
EDF.FILE.Ext = FILENAME(PPos+1:length(FILENAME));
EDF.FILE.Name = FILENAME(SPos+1:PPos-1);
if SPos==0
    EDF.FILE.Path = pwd;
else
    EDF.FILE.Path = FILENAME(1:SPos-1);
end;
EDF.FileName = [EDF.FILE.Path SLASH EDF.FILE.Name '.' EDF.FILE.Ext];

H1=char(fread(EDF.FILE.FID,256,'char')');     %
EDF.VERSION=H1(1:8);                     % 8 Byte  Versionsnummer
%if 0 fprintf(2,'LOADEDF: WARNING  Version EDF Format %i',ver); end;
EDF.PID = deblank(H1(9:88));                  % 80 Byte local patient identification
EDF.RID = deblank(H1(89:168));                % 80 Byte local recording identification
%EDF.H.StartDate = H1(169:176);         % 8 Byte
%EDF.H.StartTime = H1(177:184);         % 8 Byte
EDF.T0=[str2num(H1(168+[7 8])) str2num(H1(168+[4 5])) str2num(H1(168+[1 2])) str2num(H1(168+[9 10])) str2num(H1(168+[12 13])) str2num(H1(168+[15 16])) ];

% Y2K compatibility until year 2090
if EDF.VERSION(1)=='0'
        if EDF.T0(1) < 91
                EDF.T0(1)=2000+EDF.T0(1);
        else
                EDF.T0(1)=1900+EDF.T0(1);
        end;
else ;
        % in a future version, this is hopefully not needed
end;

EDF.HeadLen = str2num(H1(185:192));  % 8 Byte  Length of Header
% reserved = H1(193:236);            % 44 Byte
EDF.NRec = str2num(H1(237:244));     % 8 Byte  # of data records
EDF.Dur = str2num(H1(245:252));      % 8 Byte  # duration of data record in sec
EDF.NS = str2num(H1(253:256));       % 8 Byte  # of signals

EDF.Label = char(fread(EDF.FILE.FID,[16,EDF.NS],'char')');
EDF.Transducer = char(fread(EDF.FILE.FID,[80,EDF.NS],'char')');
EDF.PhysDim = char(fread(EDF.FILE.FID,[8,EDF.NS],'char')');

EDF.PhysMin= str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));
EDF.PhysMax= str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));
EDF.DigMin = str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')')); %
EDF.DigMax = str2num(char(fread(EDF.FILE.FID,[8,EDF.NS],'char')')); %

% check validity of DigMin and DigMax
if (length(EDF.DigMin) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Digital Minimum\n');
        EDF.DigMin = -(2^15)*ones(EDF.NS,1);
end
if (length(EDF.DigMax) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Digital Maximum\n');
        EDF.DigMax = (2^15-1)*ones(EDF.NS,1);
end
if (any(EDF.DigMin >= EDF.DigMax))
        fprintf(2,'Warning OPENEDF: Digital Minimum larger than Maximum\n');
end
% check validity of PhysMin and PhysMax
if (length(EDF.PhysMin) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Physical Minimum\n');
        EDF.PhysMin = EDF.DigMin;
end
if (length(EDF.PhysMax) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Physical Maximum\n');
        EDF.PhysMax = EDF.DigMax;
end
if (any(EDF.PhysMin >= EDF.PhysMax))
        fprintf(2,'Warning OPENEDF: Physical Minimum larger than Maximum\n');
        EDF.PhysMin = EDF.DigMin;
        EDF.PhysMax = EDF.DigMax;
end
EDF.PreFilt= char(fread(EDF.FILE.FID,[80,EDF.NS],'char')'); %
tmp = fread(EDF.FILE.FID,[8,EDF.NS],'char')';   %   samples per data record
EDF.SPR = str2num(char(tmp));   %   samples per data record

fseek(EDF.FILE.FID,32*EDF.NS,0);

EDF.Cal = (EDF.PhysMax-EDF.PhysMin)./ ...
    (EDF.DigMax-EDF.DigMin);
EDF.Off = EDF.PhysMin - EDF.Cal .* EDF.DigMin;
tmp = find(EDF.Cal < 0);
EDF.Cal(tmp) = ones(size(tmp));
EDF.Off(tmp) = zeros(size(tmp));

EDF.Calib=[EDF.Off';(diag(EDF.Cal))];
%EDF.Calib=sparse(diag([1; EDF.Cal]));
%EDF.Calib(1,2:EDF.NS+1)=EDF.Off';

EDF.SampleRate = EDF.SPR / EDF.Dur;

EDF.FILE.POS = ftell(EDF.FILE.FID);
if EDF.NRec == -1   % unknown record size, determine correct NRec
  fseek(EDF.FILE.FID, 0, 'eof');
  endpos = ftell(EDF.FILE.FID);
  EDF.NRec = floor((endpos - EDF.FILE.POS) / (sum(EDF.SPR) * 2));
  fseek(EDF.FILE.FID, EDF.FILE.POS, 'bof');
  H1(237:244)=sprintf('%-8i',EDF.NRec);      % write number of records
end;

EDF.Chan_Select=(EDF.SPR==max(EDF.SPR));
for k=1:EDF.NS
    if EDF.Chan_Select(k)
        EDF.ChanTyp(k)='N';
    else
        EDF.ChanTyp(k)=' ';
    end;
    if contains(upper(EDF.Label(k,:)),'ECG')
        EDF.ChanTyp(k)='C';
    elseif contains(upper(EDF.Label(k,:)),'EKG')
        EDF.ChanTyp(k)='C';
    elseif contains(upper(EDF.Label(k,:)),'EEG')
        EDF.ChanTyp(k)='E';
    elseif contains(upper(EDF.Label(k,:)),'EOG')
        EDF.ChanTyp(k)='O';
    elseif contains(upper(EDF.Label(k,:)),'EMG')
        EDF.ChanTyp(k)='M';
    end;
end;

EDF.AS.spb = sum(EDF.SPR);  % Samples per Block
bi=[0;cumsum(EDF.SPR)];

idx=[];idx2=[];
for k=1:EDF.NS,
    idx2=[idx2, (k-1)*max(EDF.SPR)+(1:EDF.SPR(k))];
end;
maxspr=max(EDF.SPR);
idx3=zeros(EDF.NS*maxspr,1);
for k=1:EDF.NS, idx3(maxspr*(k-1)+(1:maxspr))=bi(k)+ceil((1:maxspr)'/maxspr*EDF.SPR(k));end;

%EDF.AS.bi=bi;
EDF.AS.IDX2=idx2;
%EDF.AS.IDX3=idx3;


DAT.Head=EDF;
DAT.MX.ReRef=1;

%DAT.MX=feval('loadxcm',EDF);

return;
