function event = read_micromed_event(trcfile)

% reads the events of the Micromed TRC format files

fid=fopen(trcfile,'r');

%------------------reading patient & recording info----------
fseek(fid,64,-1);
surname=char(fread(fid,22,'char'))';
name=char(fread(fid,20,'char'))';

fseek(fid,128,-1);
day=fread(fid,1,'char');
if length(num2str(day))<2
    day=['0' num2str(day)];
else
    day=num2str(day);
end
month=fread(fid,1,'char');
switch month
case 1 
    month='JAN';
case 2 
    month='FEB';
case 3 
    month='MAR';
case 4 
    month='APR';
case 5 
    month='MAY';
case 6 
    month='JUN';
case 7 
    month='JUL';
case 8 
    month='AUG';
case 9 
    month='SEP';
case 10 
    month='OCT';
case 11 
    month='NOV';
case 12 
    month='DEC';
end
year=num2str(fread(fid,1,'char')+1900);

%------------------ Reading Header Info ---------

fseek(fid,175,-1);
Header_Type=fread(fid,1,'char');
if Header_Type ~= 4
    error('*.trc file is not Micromed System98 Header type 4')
end

fseek(fid,138,-1);
Data_Start_Offset=fread(fid,1,'uint32');
Num_Chan=fread(fid,1,'uint16');
Multiplexer=fread(fid,1,'uint16');
Rate_Min=fread(fid,1,'uint16');
Bytes=fread(fid,1,'uint16');
fseek(fid,176+8,-1);
Code_Area=fread(fid,1,'uint32');
Code_Area_Length=fread(fid,1,'uint32');
fseek(fid,192+8,-1);
Electrode_Area=fread(fid,1,'uint32');
Electrode_Area_Length=fread(fid,1,'uint32');

fseek(fid,400+8,-1);
Trigger_Area=fread(fid,1,'uint32');
Tigger_Area_Length=fread(fid,1,'uint32');

%----------------- Read Trace Data ----------

fseek(fid,Data_Start_Offset,-1);
switch Bytes
case 1    
    trace=fread(fid,'uint8');
case 2
    trace=fread(fid,'uint16');
case 4
    trace=fread(fid,'uint32');
end
m=length(trace);
if rem(m,Num_Chan)~=0
    roundata=floor(m/Num_Chan);
    trace=trace(1:roundata*Num_Chan);
    m=length(trace);
end
trace=reshape(trace,Num_Chan,m/Num_Chan);

%---------------- Reading Trigger Data ----------
fseek(fid,Trigger_Area,-1);
for l=1:Tigger_Area_Length/6
    trigger(1,l)=fread(fid,1,'uint32');
    trigger(2,l)=fread(fid,1,'uint16');
end

first_trigger=trigger(1,1);
m=length(trace);
tl=length(trigger);
NoTrig=0;
for tr=1:tl
    if ((trigger(1,tr) <= m) & (trigger(1,tr) >= first_trigger))
        NoTrig=NoTrig+1;
    end
end
if NoTrig > 0
   	trigger=trigger(:,1:NoTrig);
else
	trigger=[];
	first_trigger=[];
end

fclose(fid);

event = [];

if ~isempty(trigger)
  for E=1:length(trigger)
    event(E).type    = 'MARKER';
    event(E).sample  = trigger(1,E)+1;
    event(E).value   = trigger(2,E);
  end
end
