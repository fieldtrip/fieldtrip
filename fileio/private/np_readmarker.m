function [np_marker] = np_readmarker (filename, idx_begin, data_length, option)
%
% function [np_marker] = np_readmarker (filename, idx_begin, data_length, option)
%
% np_readmarker reads marker from a NEURO PRAX marker file (*.EE_).
%
% Syntax:
%
%   [np_marker] = np_readdata(filename,idx_begin,data_length,'samples');
%   [np_marker] = np_readdata(filename,idx_begin,data_length,'time');
%
% Input data:
%
%   filename    -   the complete filename with path
%                   (e. g.  C:\Document...\20030716103637.EEG)
%   idx_begin   -   the start index of the data block to be read
%   data_length -   the length of the data block to be read
%   option      -   if option = 'samples':
%                       marker will be read from sample index 'idx_begin'
%                       to sample index 'idx_begin' + 'data_length' - 1
%                   if option = 'time':
%                       marker will be read from time index 'idx_begin'
%                       to time index 'idx_begin' + 'data_length'
%
%                   To read all markers use: idx_start = 0, data_length = inf, option =
%                   'samples'.
%
% Output data:
%
%   np_marker        -   structure
%
%      np_marker.markernames   -   cell array with markernames 
%      np_marker.markertyp     -   vector array with markertypes
%      np_marker.marker        -   cell array with marker vectors 
%                                  ( = sample indices if option = 'samples',
%                                    = time indices if option = 'time');
%
% Version:  1.2. (2005-01-19)
%           1.1. (2004-10-22)
%
% 1. Artefact trials will not be considered.
% 2. Trials within pause intervals will not be considered.
%
% See also: np_readfileinfo, np_readdata
%
% eldith GmbH
% Gustav-Kirchhoff-Str. 5
% D-98693 Ilmenau
% Germany
% 22.10.2004


% Dateiinfo's lesen
np_info=np_readfileinfo(filename,'NO_MINMAX');

% Daten initialisieren, Start- und Stopsamples ermitteln
if (0==idx_begin) && (inf==data_length)
    N_start=0;
    N_stop = np_info.N - 1;                 
else
    switch upper(option)
        case 'SAMPLES', 
            N_start = idx_begin;
            N_stop  = N_start + data_length - 1;
        case 'TIME',
            N_start = round(idx_begin*np_info.fa);
            N_stop  = N_start + round(data_length*np_info.fa) - 1;
        otherwise,
            error ('Bad specification for ''option''');
    end
    if (N_start<0)
        error ('idx_begin is to small.');
    end
    if (N_stop>(np_info.N-1))
        error('data_length is to big.');
    end
end

% Markeranzahl auslesen, Strukturen und cell-Arrays definieren 
fid=fopen([filename(1:length(filename)-1) '_'],'r');
if fid==-1,
    error(['Unable to open file "' [filename(1:length(filename)-1) '_'] '" . Error code: ' ferror(fid)]);
end
try
    s=fscanf(fid,'%c',inf);
    fclose(fid);
catch
    fclose(fid);
    error('Error while reading *.EE_ file.');
end
Names=strfind(s,'Name');
Anzahl_Marker=length(Names);
np_marker=struct('marker','','markernames','');
np_marker.marker=cell(Anzahl_Marker,1);
np_marker.markertyp=zeros(1,Anzahl_Marker);
np_marker.markernames=cell(Anzahl_Marker,1);

% Markernames bestimmen und in Struktur eintragen
for i=1:Anzahl_Marker
    Idx=Names(i);                           % auf SampleIndex von "Name..." setzen
    while strcmp(s(Idx),'=')~=1,            % lesen bis Zeichen "="
        Idx=Idx+1;    
    end
    np_marker.markertyp(i)=str2num(s(Names(i)+4:Idx-1));
    Idx=Idx+1;                              
    np_marker.markernames{i}=[];            % lesen ab > "=" bis Zeilenende
    while (strcmp(s(Idx),char(13))~=1),
        np_marker.markernames{i}=[np_marker.markernames{i} s(Idx)];     % Markernamen eintragen
        Idx=Idx+1;    
    end
end

% Indexzeilen lesen und aneinanderhängen
% z.B. Index1=16518:0(Neu;Fp1#Fp1-REF_EEG;Fp2#F.....
% z.B. Index2=70:0(Neu;FB)|4:612|9:1616(100%)|3:2292|8:...
% Blöcke mit "|" einklammern
Index_Index=strfind(s,'Index');
IndexZeile=[];
for i=1:length(Index_Index)
    z=Index_Index(i);
    str=[];
    while strcmp(s(z),char(13))~=1,
        str=[str s(z)];
        z=z+1;
        if z>length(s),
            break;
        end
    end
    z=strfind(str,'=');
    IndexZeile=[IndexZeile '|' str(z+1:length(str))];
end
IndexZeile=[IndexZeile '|'];

BlockIndex=strfind(IndexZeile,'|');
Anzahl_Bloecke=length(BlockIndex)-1;

last_MarkerTyp=[];
last_MarkerName='';
last_SampleIndex=[];
last_FBMarkerTyp=[];
last_FBMarkerName='';
last_FBSampleIndex=[];

for b=1:Anzahl_Bloecke
    
    % Block extrahieren
    Block=IndexZeile(BlockIndex(b)+1:BlockIndex(b+1)-1);
   
    % MarkerTyp ermitteln (z.B. 3, 70, 16518, ...)
    % mehrer Doppelpunkte können im Block vorkommen
    doppelpunkt=find(Block==':');
    if isempty(doppelpunkt),
        MarkerTyp=last_MarkerTyp;
        doppelpunkt=0;
    else
        MarkerTyp=str2num(Block(1:doppelpunkt(1)-1));
    end
    
    % prüfen, ob es ein Artefakttrial ist
    if strcmp(Block(doppelpunkt(1)+1),'<')==1,
        Artefakttrial=1;
    else
        Artefakttrial=0;
    end
    
    % SampleIndex ermitteln, dabei nur Ziffern von 0...9
    % berücksichtigen
    %   double('0') = 48
    %   double('1') = 49
    %   ...
    %   double('9') = 57
    if Artefakttrial
        z=doppelpunkt(1)+2;
    else
        z=doppelpunkt(1)+1;
    end
    SampleIndexString=[];
    d=double(Block(z));
    while (d>=48) && (d<=57)
        SampleIndexString=[SampleIndexString Block(z)];
        z=z+1;
        if z>length(Block),
            break;
        end
        d=double(Block(z));
    end
    SampleIndex=str2num(SampleIndexString);
    
    % Behandlung des Markertyps hängt vom Markertyp ab
    % a) unbekannter Markertyp: keine Behandlung
    % b) bekannter Markertyp und MarkerTyp nicht in Menge: [1 2 3 4 6 7 8 9]
    %    --> SampleIndex in entsprechenden Markertypen eintragen
    % c) Markertyp = [6 oder 7 oder 8 oder 9] (='FBQuote+', 'FBQuote-','TRQuote+' oder 'TRQuote-')
    %    --> wenn kein Artefakt und vorhergehender Marker <> 'Pause'
    %        dann Sampleindizes für aktuellen und vorhergehenden Marker in np_marker
    %        eintragen
    %        d.h. Artefakttrials werden nicht ausgelesen, 
    %        d.h. Trials mit Pausen werden nicht gelesen
    %        d.h. Trials mit vorzeitiger Beendigung werden nicht
    %        berücksichtigt
    % Zusätzlich: letzten Feedbackmarker merken, wenn FB+, FB-, TR+ oder
    % TR-Trial gefunden wurde
    %
    % frühere Version (1.1): hier wurde der Vergleich mit:
    %       strcmp(upper(MarkerName),'FBQUOTE+')==1 oder ...
    %       durchgeführt; jetzt traten EE_-Dateien auf, in 
    %       denen anstatt "Quote" das Wort "Qoute" stand -> Fehlermeldung
    %       bzw. falsche Marker werden ausgelesen
    %       außerdem werden in Zukunft die Trialanfang- und -endemarker
    %       anders heißen, nur die Zuordnung des Markertyps bleibt;
    %       deshalb ab Version 1.2. (19.01.2005): Vergleiche über
    %       Markertypen 1-4 und 6-9 (bleiben fest).
    %       MarkerTyp 1:    TR+
    %       MarkerTyp 2:    TR-
    %       MarkerTyp 3:    FB+
    %       MarkerTyp 4:    FB-
    %       MarkerTyp 6:    TRQuote+
    %       MarkerTyp 7:    TRQuote-
    %       MarkerTyp 8:    FBQuote+
    %       MarkerTyp 9:    FBQuote-
    
    %
    % MarkerTyp kann [] sein, z. B. bei Markern, die nicht mit Namen
    % aufgelistet wurden (z. B. 16518:...); dann nächsten Block lesen
    if isempty(MarkerTyp),
        continue;
    end
    idx=find(np_marker.markertyp==MarkerTyp);
    %
    % Ergänzung: 16.2.2005
    %
    if length(idx)>1,
        error(['Typ ' num2str(MarkerTyp) ' was found ' num2str(length(idx)) ' times in np_marker.markernames.']);
    end
    if ~isempty(idx),       % Fall a) ausschließen
        MarkerName=np_marker.markernames{idx};
        if (0==ismember(MarkerTyp,[1 2 3 4 6 7 8 9])),
            % Fall b)
            if (SampleIndex>=N_start) && (SampleIndex<=N_stop)
                np_marker.marker{idx}=[np_marker.marker{idx} SampleIndex];
            end
        end

        if (1==ismember(MarkerTyp,[1 2 3 4])),
            % letzten FB-Marker merken
            last_FBMarkerTyp=MarkerTyp;
            last_FBSampleIndex=SampleIndex;
            last_FBMarkerName=MarkerName;
        end
        
        if (1==ismember(MarkerTyp,[6 7 8 9])),
            % Fall c)
            if ~Artefakttrial,
                if (strcmp(upper(MarkerName),'Pause')~=1) && (strcmp(upper(MarkerName),'PauseEnd')~=1),
                    if ~isempty(last_FBSampleIndex) && (last_FBSampleIndex>=N_start) && (last_FBSampleIndex<=N_stop)
                        % vorhergehenden Trial abspeichern
                        idx2=find(np_marker.markertyp==last_FBMarkerTyp);
                        if ~isempty(idx2),
                            np_marker.marker{idx2}=[np_marker.marker{idx2} last_FBSampleIndex];
                        end
                    end
                    if (SampleIndex>=N_start) && (SampleIndex<=N_stop)
                        % aktuellen Trial abspeichern
                        np_marker.marker{idx}=[np_marker.marker{idx} SampleIndex];
                    end % if SampleIndex>=Nstart
                end % if ~Pause
            end % if ~Artefakttrial
        end % if Markername = Endmarker FB/TR
        
        % aktualisiere letzten Marker
        last_MarkerName=MarkerName;
        last_MarkerTyp=MarkerTyp;
        last_SampleIndex=SampleIndex;
    end % if ~isempty...
end % for b=1:...
