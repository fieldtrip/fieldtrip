function [HDR]=openxml(arg1,CHAN,arg4,arg5,arg6)
% OPENXML reads XML files and tries to extract biosignal data
%
% This is an auxilary function to SOPEN. 
% Use SOPEN instead of OPENXML.
%

%
% HDR = openxml(HDR);
%
% HDR contains the Headerinformation and internal data
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	Copyright 2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

if ischar(arg1); 
        HDR.FileName = arg1; 
        HDR.FILE.PERMISSION = 'r'; 
else
        HDR = arg1; 
end;         

%if strncmp(HDR.TYPE,'XML',3),
%        if any(HDR.FILE.PERMISSION=='r'),
                fid = fopen(HDR.FileName,HDR.FILE.PERMISSION,'ieee-le');
                s = char(fread(fid,[1,1024],'char')); 
                if all(s(1:2)==[255,254]) & all(s(4:2:end)==0)
                        HDR.TYPE='XML-UTF16';
                elseif ~isempty(findstr(char(s),'?xml version'))
                        HDR.TYPE='XML-UTF8';
		end; 
		fseek(fid,0,'bof');
                if strcmp(HDR.TYPE,'XML-UTF16'),
                        magic = char(fread(fid,1,'uint16'));
                        HDR.XML = char(fread(fid,[1,inf],'uint16'));
                elseif strcmp(HDR.TYPE,'XML-UTF8'),
                        HDR.XML = char(fread(fid,[1,inf],'uint8'));
                end;
                fclose(fid);
                HDR.FILE.FID = fid;
                

                if 0, ~exist('xmlstruct')
                        warning('XML toolbox missing')
                end;
                if 1, 
                        HDR.XMLstruct = xmlstruct(HDR.XML,'sub');
                        HDR.XMLlist   = xmlstruct(HDR.XML);
                end;

                try
                        XML = xmltree(HDR.XML);
                        XML = convert(XML);
                        HDR.XML  =  XML; 
			HDR.TYPE = 'XML';
                catch
                        fprintf(HDR.FILE.stderr,'ERROR SOPEN (XML): XML-toolbox missing or invalid XML file.\n');
                        return;
                end;

                tmp = fieldnames(HDR.XML);
                if any(strmatch('PatientDemographics',tmp)) & any(strmatch('TestDemographics',tmp)) & any(strmatch('RestingECGMeasurements',tmp)) & any(strmatch('Diagnosis',tmp)) & any(strmatch('Waveform',tmp)),
                        % GE-Marquette FDA-XML   MAC5000 
                        
                        HDR.Patient.ID = HDR.XML.PatientDemographics.PatientID;
                        tmp = HDR.XML.PatientDemographics.Gender;
                        HDR.Patient.Sex = 2*(upper(tmp(1))=='F') + (upper(tmp(1))=='M');
                        HDR.Patient.Name = [HDR.XML.PatientDemographics.PatientLastName,', ',HDR.XML.PatientDemographics.PatientFirstName];

                        tmp = HDR.XML.TestDemographics.AcquisitionDate; 
                        tmp(tmp=='-') = ' '; 
                        HDR.T0([3,2,1])=str2double(tmp);
                        tmp = HDR.XML.TestDemographics.AcquisitionTime; 
                        tmp(tmp==':') = ' '; 
                        HDR.T0(4:6) = str2double(tmp);

                        HDR.NS = str2double(HDR.XML.Waveform.NumberofLeads);
                        HDR.SampleRate = str2double(HDR.XML.Waveform.SampleBase); 
                        HDR.Filter.LowPass  = str2double(HDR.XML.Waveform.LowPassFilter)*ones(1,HDR.NS);
                        HDR.Filter.HighPass = str2double(HDR.XML.Waveform.HighPassFilter)*ones(1,HDR.NS);
                        HDR.Filter.Notch    = str2double(HDR.XML.Waveform.ACFilter)*ones(1,HDR.NS);

                        HDR.NRec = 1; 
                        HDR.SPR = 1;
                        for k = 1:HDR.NS,
                                CH = HDR.XML.Waveform.LeadData{k};
                                HDR.AS.SPR(k)  = str2double(CH.LeadSampleCountTotal);
                                HDR.SPR        = lcm(HDR.SPR,HDR.AS.SPR(k));
                                HDR.Cal(k)     = str2double(CH.LeadAmplitudeUnitsPerBit);
                                HDR.PhysDim{k} = CH.LeadAmplitudeUnits;
                                HDR.Label{k}   = CH.LeadID;
                        
                                t = radix64d(CH.WaveFormData); 
                                t = 256*t(2:2:end) + t(1:2:end);
                                t = t - (t>=(2^15))*(2^16);
                                HDR.data(:,k) = t(:); 
                        end;
                        HDR.TYPE = 'native'; 

                        
                elseif any(strmatch('component',tmp))
                          % FDA-XML Format
                        tmp   = HDR.XML.component.series.derivation;
                        if isfield(tmp,'Series');
                                tmp = tmp.Series.component.sequenceSet.component;
                        else    % Dovermed.CO.IL version of format
                                tmp = tmp.derivedSeries.component.sequenceSet.component;
                        end;
                        HDR.NS = length(tmp)-1;
                        HDR.NRec = 1; 
                        HDR.Cal = 1;
                        HDR.PhysDim = {' '};
                        HDR.SampleRate = 1;
                        HDR.TYPE = 'XML-FDA';     % that's an FDA XML file 

                        
                elseif any(strmatch('dataacquisition',tmp)) & any(strmatch('reportinfo',tmp)) & any(strmatch('patient',tmp)) & any(strmatch('documentinfo',tmp)),
                        % SierraECG  1.03  *.open.xml from PHILIPS
                        HDR.SampleRate = str2double(HDR.XML.dataacquisition.signalcharacteristics.samplingrate);
                        HDR.NS  = str2double(HDR.XML.dataacquisition.signalcharacteristics.numberchannelsvalid);
                        HDR.Cal = str2double(HDR.XML.reportinfo.reportgain.amplitudegain.overallgain);
                        HDR.PhysDim = {'uV'};
                        HDR.Filter.HighPass = str2double(HDR.XML.reportinfo.reportbandwidth.highpassfiltersetting);
                        HDR.Filter.LowPass  = str2double(HDR.XML.reportinfo.reportbandwidth.lowpassfiltersetting);
                        HDR.Filter.Notch    = str2double(HDR.XML.reportinfo.reportbandwidth.notchfiltersetting);

                        t = HDR.XML.reportinfo.reportformat.waveformformat.mainwaveformformat;
                        k = 0;
                        HDR.Label={};
                        while ~isempty(t),
                                [s,t] = strtok(t,' ');
                                k = k+1;
                                HDR.Label{k, 1} = [s,' '];
                        end;
                        HDR.Patient.Id     = str2double(HDR.XML.patient.generalpatientdata.patientid);
                        tmp    = HDR.XML.patient.generalpatientdata.age;
                        if isfield(tmp,'years'),
                                HDR.Patient.Age    = str2double(tmp.years);
                        end
                        if isfield(tmp,'dateofbirth')
                                tmp = tmp.dateofbirth;
                                tmp(tmp=='-')=' ';
                                HDR.Patient.Birthday([6,5,4]) = str2double(tmp);
                        end;

                        tmp    = HDR.XML.patient.generalpatientdata.sex;
                        HDR.Patient.Sex    = strncmpi(tmp,'Male',1) + strncmpi(tmp,'Female',1)*2;
                        HDR.Patient.Weight = str2double(HDR.XML.patient.generalpatientdata.weight.kg);
                        HDR.Patient.Height = str2double(HDR.XML.patient.generalpatientdata.height.cm);

                        HDR.VERSION = HDR.XML.documentinfo.documentversion;
                        HDR.TYPE = HDR.XML.documentinfo.documenttype;


                elseif any(strmatch('component',tmp)) & any(strmatch('reportinfo',tmp)) & any(strmatch('patient',tmp)) & any(strmatch('documentinfo',tmp)),
                           % FDA-XML Format
                        tmp   = HDR.XML.component.series.derivation;
                        if isfield(tmp,'Series');
                                tmp = tmp.Series.component.sequenceSet.component;
                        else    % Dovermed.CO.IL version of format
                                tmp = tmp.derivedSeries.component.sequenceSet.component;
                        end;
                        HDR.NS = length(tmp)-1;
                        HDR.NRec = 1; 
                        HDR.Cal = 1;
                        HDR.PhysDim = {' '};
                        HDR.SampleRate = 1;
                        HDR.TYPE = 'XML-FDA';     % that's an FDA XML file 


                elseif any(strmatch('StripData',tmp)) & any(strmatch('ClinicalInfo',tmp)) & any(strmatch('PatientInfo',tmp)) & any(strmatch('ArrhythmiaData',tmp)),
                        % GE Case8000 stress ECG

                	HDR.SampleRate = str2double(HDR.XML.StripData.SampleRate);
                	tmp = HDR.XML.ClinicalInfo.ObservationDateTime; 
                	HDR.T0 = [str2double(tmp.Year), str2double(tmp.Month), str2double(tmp.Day), str2double(tmp.Hour), str2double(tmp.Minute), str2double(tmp.Second)]; 

                	HDR.Patient.Id = HDR.XML.PatientInfo.PID; 
                	HDR.Patient.Name = 'X'; % [HDR.XML.PatientInfo.Name, ', ', HDR.XML.PatientInfo.GivenName];
                	HDR.Patient.Age = str2double(HDR.XML.PatientInfo.Age);
                	tmp = HDR.XML.PatientInfo.Gender; 
                	HDR.Patient.Sex = any(tmp(1)=='Mm') + any(tmp(1)=='Ff')*2;
                	HDR.Patient.Height = str2double(HDR.XML.PatientInfo.Height);
                	HDR.Patient.Weight = str2double(HDR.XML.PatientInfo.Weight); 
                	tmp = HDR.XML.PatientInfo.BirthDateTime; 
                	HDR.Patient.Birthday = [str2double(tmp.Year), str2double(tmp.Month), str2double(tmp.Day),0,0,0]; 
                	
                	tmp = HDR.XML.StripData.Strip; 
                	HDR.NS = length(tmp{1}.WaveformData);
               		tmax   = str2double(tmp{end}.Time.Minute)*60 + str2double(tmp{end}.Time.Second)+10;
                	data   = repmat(NaN,tmax*HDR.SampleRate,HDR.NS); 
                	for k  = 1:length(tmp);
	                	t = HDR.SampleRate*(str2double(tmp{k}.Time.Minute)*60 + str2double(tmp{k}.Time.Second));
	                	for k2 = 1:HDR.NS, 
		                	x = str2double(tmp{k}.WaveformData{k2}); 
		                	data(t+1:t+length(x),k2)=x(:); 
	                	end; 
                	end;
                	tmp = HDR.XML.ArrhythmiaData.Strip; 
                	for k  = 1:length(tmp);
	                	t = HDR.SampleRate*(str2double(tmp{k}.Time.Minute)*60 + str2double(tmp{k}.Time.Second));
	                	for k2 = 1:HDR.NS, 
		                	x = str2double(tmp{k}.WaveformData{k2}); 
		                	data(t+1:t+length(x),k2)=x(:); 
	                	end; 
                	end;
                	HDR.data = data - 2^12*(data>2^11); 
                	HDR.TYPE = 'native'; 
                	HDR.NRec = 1; 
                	HDR.SPR  = size(HDR.data,1); 
                	HDR.Calib = sparse(2:HDR.NS,1:HDR.NS,1); 
                	HDR.FLAG.UCAL = 1; 

                else
                        fprintf(HDR.FILE.stderr,'Warning SOPEN (XML): File %s is not supported.\n',HDR.FileName);
                        return;
                end

                
                try
                	tmp=HDR.XML.componentOf.timepointEvent.componentOf.subjectAssignment.subject.trialSubject.subjectDemographicPerson.name; 
                	HDR.Patient.Name = sprintf('%s, %s',tmp.family, tmp.given);
		catch
                end;        
                
                
                HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,HDR.Cal);
                HDR.FILE.OPEN = 1;
                HDR.FILE.POS  = 0;
%        end;
%end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Auxillary functions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = radix64d(x);
% RADIX64D - decoding of radix64 encoded sequence


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%	$Id$
%	(C) 2006 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


global BIOSIG_GLOBAL

if ~isfield(BIOSIG_GLOBAL,'R64E');
	BIOSIG_GLOBAL.R64E = ['A':'Z','a':'z','0':'9','+','/'];
	BIOSIG_GLOBAL.R64D = zeros(256,1)-1;
	for k = 1:length(BIOSIG_GLOBAL.R64E),
		BIOSIG_GLOBAL.R64D(BIOSIG_GLOBAL.R64E(k)) = k-1; 
	end;
end; 

% http://www.faqs.org/rfcs/rfc2440.html

t = BIOSIG_GLOBAL.R64D(x);
t = [t(t>=0); zeros(3,1)];
N = floor(length(t)/4);
t = reshape(bitand(t(1:N*4),2^6-1),4,N);

y(1,:) = bitshift(t(1,:),2)         + bitshift(t(2,:),-4);
y(2,:) = bitshift(mod(t(2,:),16),4) + bitshift(t(3,:),-2); 
y(3,:) = bitshift(mod(t(3,:),4),6)  + t(4,:);
 
y = y(:)';
y = y(1:end-sum(x=='=')); % remove possible pad characters

