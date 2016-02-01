function [np_info] = np_read_splitted_fileinfo (filename, option)
%
% function [np_info] = np_read_splitted_fileinfo (filename, option)
%
% This function is necessary for np_readfileinfo.m.
%
% eldith GmbH
% Gustav-Kirchhoff-Str. 5
% D-98693 Ilmenau
% Germany
% 02.02.2005



% -------------------------------------------------------------------------
% Initialisierung der Struktur
% -------------------------------------------------------------------------
np_info=struct('filename','','pathname','','name','','firstname','','birthday','','ID',0,'date','00.00.00',...
               'time','00:00:00','duration','','setup','','pmtype','','algorithm','','channels','','channeltypes','',...
               'units','','PhysMin',0,'PhysMax',0,'N',0,'K',0,'fa',0,'fp_data',0);

% -------------------------------------------------------------------------
% Auslesen: fp_text, fp_data, K, N, fa
% -------------------------------------------------------------------------
d=dir(filename);
if d(1).bytes==0,
    error('File size = 0 KB.');
end

fid=fopen(filename,'r');
if fid==-1,
    error(['Unable to open file "' filename '". Error code: ' ferror(fid)]);
end

status=fseek(fid,16,'bof');
if status~=0,
    fclose(fid);
    error('Unable to set filepointer to index 16.');
end
[np_info.fp_data,count]=fread(fid,1,'int32');
if count~=1,
    fclose(fid);
    error('Unable to read filepointer for data block.');
end
[fp_header,count]=fread(fid,1,'int32');
if count~=1,
    fclose(fid);
    error('Unable to read filepointer for header block.');
end
[fp_marker,count]=fread(fid,1,'int32');
if count~=1,
    fclose(fid);
    error('Unable to read filepointer for marker block.');
end
[fp_text,count]=fread(fid,1,'int32');
if count~=1,
    fclose(fid);
    error('Unable to read filepointer for text block.');
end
status=fseek(fid,fp_header,'bof');
if status~=0,
    fclose(fid);
    error('Unable to set file to fp_header.');
end
[np_info.K,count]=fread(fid,1,'short');                     % lese Kanalanzahl
if count~=1,
    fclose(fid);
    error('Unable to read number of channels.');
end
[ignore,count]=fscanf(fid,'%19c',1);
if count~=1,
    fclose(fid);
    error('Unable to read duration string.');
end
% Sampleanzahl ermitteln:
% frühere Version:
%       [np_info.N,count]=fread(fid,1,'int32');                
%       if count~=1,
%           error(['Fehler beim Lesen der Messdauer. Fehlercode: ' ferror(fid)]);
%       end
% jetzt: Ermittlung aus FilePointern und FileSize!!!
d=dir(filename);
filesize=d.bytes;
SIZEOF_FLOAT=4;
fp_vector=sort([fp_text fp_marker fp_header np_info.fp_data filesize]);
idx=find(fp_vector==np_info.fp_data);
%
% Ergänzung: 16.2.05
% wenn fp_data = filesize, dann hat die Datei 
% wahrscheinlich keine Daten
%
if (np_info.fp_data==filesize),
    error(['Unable to read data (no data?). fp_text = ' num2str(fp_text) '; fp_marker = ' num2str(fp_marker) '; fp_header = ' num2str(fp_header) ';fp_data = ' num2str(np_info.fp_data) '; filesize = ' num2str(filesize)]);
end
np_info.N=(fp_vector(idx+1)-fp_vector(idx))/SIZEOF_FLOAT/np_info.K;
status=fseek(fid,fp_header+39,'bof');
if status~=0,
    fclose(fid);
    error('Unable to set filepointer to sampling frquency of 1st channel.');
end
[fa_str,count]=fscanf(fid,'%19c',1);                    % lese Abtastfrequenz Kanal 1
if count~=1,
    fclose(fid);
    error('Unable to read sampling frquency of 1st channel.');
end
np_info.fa=str2num(fa_str(2:19));                           % Umwandlung: CHAR-Array in double-Zahl

% -------------------------------------------------------------------------
% Auslesen: Name, Vorname, Geburtstag, ID
% -------------------------------------------------------------------------
if fp_text>0,
    status=fseek(fid,fp_text,'bof');
    if status~=0,
        fclose(fid);
        error('Unable to set filepointer to begin of text block.');
    end
    tline=fgetl(fid);       %  [PatInfo]
    tline=fgetl(fid);       %  Name=...
    np_info.name=tline(6:length(tline));
    tline=fgetl(fid);       %  Vorname=...
    np_info.firstname=tline(9:length(tline));
    tline=fgetl(fid);       %  GebDat=...
     da=tline(8:length(tline));
    np_info.birthday=da;
    %
    % alte Version hatte nur Probleme, wenn unterschiedliche Regions-
    % und Sprachoptionen -> deshalb String für GebDat nicht formatieren
    % 07.04.2005
%     if ~isempty(strfind(da,'/')),            % amerikanisches Datumsformat?
%         %bc=datevec(da);                      % = frühere Version   
%         bc=datevec(da,'dd/mm/yyyy');
%         bs=datestr(bc,24);
%         np_info.birthday=regexprep(bs,'/','.');
%     else
%         np_info.birthday=da;     % sonst Annahme deutsches/europäisches Format
%     end
    tline=fgetl(fid);       %  ID=...
    np_info.ID=tline(4:length(tline));
end

% -------------------------------------------------------------------------
% Kanalbezeichnungen lesen
% -------------------------------------------------------------------------
np_info.channels=cell(1,np_info.K);
for i=1:np_info.K
    status=fseek(fid,fp_header+35+(i-1)*203+166,'bof');
    if status~=0,
        fclose(fid);
        error(['Unable to set filepointer to signal name of channel ' num2str(i) '.']);
    end
    [h,count]=fscanf(fid,'%8c',1);
    if count~=1,
        fclose(fid);
        error(['Unable to read channel name for channel ' num2str(i) '.']);
    end
    %
    % wahre Stringlänge bestimmen!
    % Änderung am 09.12.2004
    %
       hh=h(2:8);
       laenge=0;
       for k=1:length(hh)
          if (1==strcmp(hh(k),char(0))),
               break;
          else
               laenge=laenge+1;
          end
       end
    %
    np_info.channels{i}=hh(1:laenge);
end

% -------------------------------------------------------------------------
% Kanaltypen lesen
% -------------------------------------------------------------------------
np_info.channeltypes=cell(1,np_info.K);
for i=1:np_info.K
    status=fseek(fid,fp_header+35+(i-1)*203+186,'bof');
    if status~=0,
        fclose(fid);
        error(['Unable to set filepointer to unit of channel ' num2str(i) '.']);
    end
    [h,count]=fscanf(fid,'%16c',1);
    if count~=1,
        fclose(fid);
        error(['unable to read unit of channel ' num2str(i) '.']);
    end
    h_unit=h(2:8);
    %
    % wahre Stringlängen bestimmen!
    % Änderung am 01.02.2005
    %
    laenge=0;
    for k=1:length(h_unit)
       if (1==strcmp(h_unit(k),char(0))),
            break;
       else
            laenge=laenge+1;
       end
    end
    np_info.units{i}=h_unit(1:laenge);

    h_type=h(10:16);
    %
    % wahre Stringlängen bestimmen!
    % Änderung am 01.02.2005
    %
    laenge=0;
    for k=1:length(h_type)
       if (1==strcmp(h_type(k),char(0))),
            break;
       else
            laenge=laenge+1;
       end
    end
    np_info.channeltypes{i}=h_type(1:laenge);
end
fclose(fid);

% -------------------------------------------------------------------------
% Setup, PMTyp und ALgorithmus lesen
% -------------------------------------------------------------------------
% --- alte Version: fid=fopen([filename(1:length(filename)-1) '_'],'r');
% --- bis 19.09.2005: in alter Version werden nicht die *X*.EEG Dateien
% --- für EEG-Splitting benutzt; neue Version liest Setup, PMtype,
% --- Algorithmus immer von der "ersten" EE_-Datei
[pa,fn,ex]=fileparts(filename);
if strcmp(pa,''),
    pa=pwd;
end
fid=fopen([pa filesep fn(1:14) '.EE_'],'r');
if fid==-1,
    error('Unable to read setup (marker 70) in *.EE_ file.');
end
s=fscanf(fid,'%c',inf);
fclose(fid);
Idx1=strfind(s,'70:');      % Name70=PMTyp
if ~isempty(Idx1),
    % der Marker:   70:0(...) ist vorhanden
    % Beispiel:     70:0(EEG-4-FB;FB;TP_ADHD_SCP)
    [block,R]=strtok(s(Idx1(1):length(s)),'|');    % nächsten Block bis '|' ermitteln
    [A,B]=strtok(block,'(');                       % suche nach '(' in block
    [C,D]=strtok(B(2:length(B)),')');              % suche nach ')' in Restblock C
    [np_info.setup,R]=strtok(C,';');               % Abtrennung für Setup
    [np_info.pmtype,R]=strtok(R,';');              % Abtrennung für PMTyp
    [np_info.algorithm,R]=strtok(R,';');           % Abtrennung für Algorithmus
else
    % der Marker:   70:0(...) ist nicht vorhanden
    % es muss 16518:... benutzt werden
    % nur setup kann ermittelt werden aus der letzten
    % Sekundärmontage
    fid=fopen([filename(1:length(filename)-1) '_'],'r');
    if fid==-1,
        error('Unable to read setup (marker 16518) in *.EE_ file.');
    end
    s=fscanf(fid,'%c',inf);
    fclose(fid);
    MarkerIdx=strfind(s,'16518:');       % alle Marker mit 16518: finden
    letzteSekMontage='';
    for i=1:length(MarkerIdx)
        % Beispiel 1: 16518:0|              -> keine Montage vorhanden
        % Beispiel 2: 16518:0(EEG;Fp1...)
        [block,R]=strtok(s(MarkerIdx(i):length(s)),'|');    % nächsten Block bis '|' ermitteln
        [A,B]=strtok(block,'(');                            % im Block nach '(' suchen
        [C,D]=strtok(B(2:length(B)),')');                   % im Restblock B nach ')' suchen
        [E,F]=strtok(C,';');                                % im Block C nach ';' suchen
        if length(E)>5,
            if strcmp(E(1:5),'Prim.')~=1,               % Prim.mont. oder Prim.mont.;
                letzteSekMontage=E;
            end
        else
            letzteSekMontage=E;
        end
    end   
    np_info.setup=letzteSekMontage;
end
% Ergänzung: 21.09.05
% Alias für Algorithmus (Protokoll) eintragen
dw=which('np_read_splitted_fileinfo.m');
if ~isempty(dw),
    [pa,fn,ex]=fileparts(dw);
    if strcmp(pa,''),
        pa=pwd;
    end
    try
        algo=load([pa filesep 'npalgo.mat']);
        idx=find(strcmp(algo.algo(:,1),np_info.algorithm)==1);
        if ~isempty(idx),
            np_info.algorithm=cell2mat(algo.algo(idx,2));
        end
    catch
    end
end

% -------------------------------------------------------------------------
% Dateinamen, Datum, Uhrzeit und Messdauer speichern
% -------------------------------------------------------------------------
[np_info.pathname,fn,ext]=fileparts(filename);
if strcmp(np_info.pathname,'')
    np_info.pathname=pwd;
end
np_info.filename=[fn ext];
[pa,fn,ex]=fileparts(filename);
s=fn(1:14);     % nur ein Datum für gesplittete EEG-Dateien
% 
% Anpassung an ISO 8601  'yyyy-mm-dd'             2000-03-01
% alte Version: np_info.date=datestr([s(5:6) '/' s(7:8) '/' s(1:4)]);
np_info.date=[s(1:4) '-' s(5:6) '-' s(7:8)];
np_info.duration=np_info.N/np_info.fa;
%
% neue Version: Auslesen der Zeit aus der EE_ Datei
%
fid=fopen([filename(1:length(filename)-1) '_'],'r');
if fid==-1,
    error('Unable to read primary setup in *.EE_ file.');
end
s=fscanf(fid,'%c',inf);
fclose(fid);
Idx1=strfind(s,'16384:0');
if length(Idx1)>1,      % für den Fall, dass 16384:0 mehrfach auftaucht
    Idx1=Idx1(1);       % das kann bei Pausen der Fall sein
end
if ~isempty(s)
    %
    % Formate: 16384:0(01.02.2005 12:54:41)  oder  16384:0(2/1/2005 8:50:03 PM)
    % in einer neuen Version: 01_02_2005 12:54:41
    %
    [block,R]=strtok(s(Idx1:length(s)),'|');            % nächsten Block bis '|' ermitteln
    [A,B]=strtok(block,'(');                            % im Block nach '(' suchen
    [C,D]=strtok(B(2:length(B)),')');                   % im Restblock B nach ')' suchen    
    [str,R]=strtok(C,' ');                              % ergibt das Datum in str und time in R
    np_info.time=datestr(datenum(R(2:length(R))),13);   % formatiere Zeit: hh:mm:ss
    %
    % Ergänzung 07.04.2005:
    % Prüfen des Messdatums
    % z.B. EEG-Dateiname: 20050407235920.EEG   07.04.2005 23:59:20
    %      Messbeginn:    8.4.2005 00:20:21
    % dann muss das Messdatum um 1 Tag korrigiert werden
    % Alternative: Messdatum aus EE_-Datei -> abhängig von Regions- und
    % Spracheinstellungen
    filename_time=str2num(np_info.filename(9:14));          % z.B. 235920
    info_time=str2num([np_info.time(1:2) np_info.time(4:5) np_info.time(7:8)]);     % z.B. 002021
    if (info_time<filename_time),
        % Messgebinn war nach Mitternacht, Messdatum um 1 Tag erhöhen
        np_info.date=datestr(datenum(str2num(np_info.date(1:4)),...             % Jahr
                                     str2num(np_info.date(6:7)),...             % Monat 
                                     str2num(np_info.date(9:10)))+1,29);        % Tag um 1 erhöhen
    end
end

% -------------------------------------------------------------------------
% Physikalisches Minimum und Maximum ermitteln
% -------------------------------------------------------------------------
np_info.PhysMin=zeros(1,np_info.K);
np_info.PhysMax=zeros(1,np_info.K);
if (nargin==2) && (strcmp(upper(option),'NO_MINMAX'))
    return;
end
fid=fopen([np_info.pathname filesep np_info.filename],'r');
if fid==-1,
    error('Error while opening *.EEG file (read PhysMinMax).');
end
status=fseek(fid,np_info.fp_data,'bof');
if status~=0,
    fclose(fid);
    error('Unable to set filepointer to begin of data block (read PhysMinMax).');
end
SamplesToRead=1000;
try
    for k=1:np_info.K                   % setze Werte für PhysMin und PhysMax
        val=fread(fid,1,'float');       % auf die Anfangswerte des jeweiligen
        np_info.PhysMin(1,k)=val;       % Kanals
        np_info.PhysMax(1,k)=val;
    end
    while 1,
        [val,count]=fread(fid,[np_info.K SamplesToRead],'float');
        val=val';
        np_info.PhysMin=min([min(val,[],1); np_info.PhysMin],[],1);
        np_info.PhysMax=max([max(val,[],1); np_info.PhysMax],[],1);
        if count<(np_info.K*SamplesToRead),
            break;
        end
    end
    % wichtige Ergänzung: 23.08.2004
    fclose(fid);
catch
    fclose(fid);
    error('Error while calculating physical minimum or maximum.');
end
