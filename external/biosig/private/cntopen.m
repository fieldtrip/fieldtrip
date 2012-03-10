function [CNT,h,e]=cntopen(arg1,arg2,arg3,arg4,arg5,arg6)
% CNTOPEN opens neuroscan files (but does not read the data). 
% However, it is recommended to use SOPEN instead of CNTOPEN.
% For loading whole Neuroscan data files, use SLOAD. 
%
% see also: SLOAD, SOPEN, SREAD, SCLOSE, SEOF, STELL, SSEEK.

%	$Id$
%	Copyright (c) 1997-2006,2007,2008,2009 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.

if nargin<2, arg2=''; end; 
if isstruct(arg1),
	CNT=arg1;
	if CNT.FILE.OPEN,
                status=fseek(CNT.FILE.FID,0,'bof');	
                if status,
                        fprintf(CNT.FILE.stderr,'Warning CNTOPEN: I/O error in file %s\n',CNT.FileName);
                end;        
        else
		CNT.FILE.FID = fopen(CNT.FileName,CNT.FILE.PERMISSION,'ieee-le');          
		CNT.FILE.OPEN= 1;
	end;
else
	CNT.FileName = arg1;
        [CNT.FILE.Path,CNT.FILE.Name,ext]=fileparts(arg1);
        CNT.FILE.Ext = ext(2:end);
        CNT.FILE.FID = fopen(CNT.FileName,CNT.FILE.PERMISSION,'ieee-le');          
	if CNT.FILE.FID<0,
		fprintf(2,'Error CNTOPEN: file %s not found.\n',CNT.FileName); 
		return;
	end;
	CNT.FILE.OPEN = 1;
end;

fid = CNT.FILE.FID;

%%%%% READ HEADER
if 0,   % old header
        %h.rev               = fread(fid,12,'uint8');
        %h.nextfile          = fread(fid,1,'long');
        %h.prevfile          = fread(fid,1,'long');
        h.type              = fread(fid,1,'uint8');
        h.id                = fread(fid,20,'uint8');
        h.oper              = fread(fid,20,'uint8');
        h.doctor            = fread(fid,20,'uint8');
        h.referral          = fread(fid,20,'uint8');
        h.hospital          = fread(fid,20,'uint8');
        h.patient           = fread(fid,20,'uint8');
        h.age               = fread(fid,1,'short');
        h.sex               = fread(fid,1,'uint8');
        h.hand              = fread(fid,1,'uint8');
        h.med               = fread(fid,20,'uint8');
        h.category          = fread(fid,20,'uint8');
        h.state             = fread(fid,20,'uint8');
        h.label             = fread(fid,20,'uint8');
        h.date              = fread(fid,10,'uint8');
        h.time              = fread(fid,12,'uint8');
        h.avgmode           = fread(fid,1,'uint8');
        h.review            = fread(fid,1,'uint8');
        h.nsweeps           = fread(fid,1,'ushort');
        h.compsweeps        = fread(fid,1,'ushort');
        h.pnts              = fread(fid,1,'ushort');
        h.nchannels         = fread(fid,1,'short');
        h.avgupdate         = fread(fid,1,'short');
        h.domain            = fread(fid,1,'uint8');
        h.rate              = fread(fid,1,'ushort');
        h.scale             = fread(fid,1,'double');
        h.veogcorrect       = fread(fid,1,'uint8');
        h.veogtrig          = fread(fid,1,'float');
        h.veogchnl          = fread(fid,1,'short');
        h.heogcorrect       = fread(fid,1,'uint8');
        h.heogtrig          = fread(fid,1,'float');
        h.heogchnl          = fread(fid,1,'short');
        h.baseline          = fread(fid,1,'uint8');
        h.offstart          = fread(fid,1,'float');
        h.offstop           = fread(fid,1,'float');
        h.reject            = fread(fid,1,'uint8');
        h.rejchnl1          = fread(fid,1,'uint8');
        h.rejchnl2          = fread(fid,1,'uint8');
        h.rejchnl3          = fread(fid,1,'uint8');
        h.rejchnl4          = fread(fid,1,'uint8');
        h.rejstart          = fread(fid,1,'float');
        h.rejstop           = fread(fid,1,'float');
        h.rejmin            = fread(fid,1,'float');
        h.rejmax            = fread(fid,1,'float');
        h.trigtype          = fread(fid,1,'uint8');
        h.trigval           = fread(fid,1,'float');
        h.trigchnl          = fread(fid,1,'uint8');
        h.trigisi           = fread(fid,1,'float');
        h.trigmin           = fread(fid,1,'float');
        h.trigmax           = fread(fid,1,'float');
        h.trigdur           = fread(fid,1,'float');
        h.dir               = fread(fid,1,'uint8');
        h.dispmin           = fread(fid,1,'float');
        h.dispmax           = fread(fid,1,'float');
        h.xmin              = fread(fid,1,'float');
        h.xmax              = fread(fid,1,'float');
        h.ymin              = fread(fid,1,'float');
        h.ymax              = fread(fid,1,'float');
        h.zmin              = fread(fid,1,'float');
        h.zmax              = fread(fid,1,'float');
        h.lowcut            = fread(fid,1,'float');
        h.highcut           = fread(fid,1,'float');
        h.common            = fread(fid,1,'uint8');
        h.savemode          = fread(fid,1,'uint8');
        h.manmode           = fread(fid,1,'uint8');
        h.ref               = fread(fid,20,'uint8');
        h.screen            = fread(fid,80,'uint8');
        h.seqfile           = fread(fid,80,'uint8');
        h.montage           = fread(fid,80,'uint8');
        h.heegcorrect       = fread(fid,1,'uint8');
        h.variance          = fread(fid,1,'uint8');
        h.acceptcnt         = fread(fid,1,'ushort');
        h.rejectcnt         = fread(fid,1,'ushort');
        h.reserved74        = fread(fid,74,'uint8');
       
        for n = 1:64,%h.nchannels
                e(n).lab            = fread(fid,10,'uint8');
                e(n).x_coord        = fread(fid,1,'float');
                e(n).y_coord        = fread(fid,1,'float');
                e(n).alpha_wt       = fread(fid,1,'float');
                e(n).beta_wt        = fread(fid,1,'float');
        end

        
else    % new header
        h.rev               = fread(fid,12,'uint8');
        h.nextfile          = fread(fid,1,'uint32');
        h.prevfile          = fread(fid,1,'uint32');
        h.type              = fread(fid,1,'uint8');
        h.id                = fread(fid,20,'uint8');
        h.oper              = fread(fid,20,'uint8');
        h.doctor            = fread(fid,20,'uint8');
        h.referral          = fread(fid,20,'uint8');
        h.hospital          = fread(fid,20,'uint8');
        h.patient           = fread(fid,20,'uint8');
        h.age               = fread(fid,1,'short');
        h.sex               = fread(fid,1,'uint8');
        h.hand              = fread(fid,1,'uint8');
        h.med               = fread(fid,20,'uint8');
        h.category          = fread(fid,20,'uint8');
        h.state             = fread(fid,20,'uint8');
        h.label             = fread(fid,20,'uint8');
        h.date              = fread(fid,10,'uint8');	%%%
        h.time              = fread(fid,12,'uint8');	%%%
        h.mean_age          = fread(fid,1,'float');
        h.stdev             = fread(fid,1,'float');
        h.n                 = fread(fid,1,'short');
        h.compfile          = fread(fid,38,'uint8');
        h.spectwincomp      = fread(fid,1,'float');
        h.meanaccuracy      = fread(fid,1,'float');
        h.meanlatency       = fread(fid,1,'float');
        h.sortfile          = fread(fid,46,'uint8');
        h.numevents         = fread(fid,1,'int');	%%%
        h.compoper          = fread(fid,1,'uint8');
        h.avgmode           = fread(fid,1,'uint8');
        h.review            = fread(fid,1,'uint8');
        h.nsweeps           = fread(fid,1,'ushort');
        h.compsweeps        = fread(fid,1,'ushort');
        h.acceptcnt         = fread(fid,1,'ushort');
        h.rejectcnt         = fread(fid,1,'ushort');
        h.pnts              = fread(fid,1,'ushort');
        h.nchannels         = fread(fid,1,'ushort');	%%%
        h.avgupdate         = fread(fid,1,'ushort');
        h.domain            = fread(fid,1,'uint8');
        h.variance          = fread(fid,1,'uint8');
        h.rate              = fread(fid,1,'ushort');	%%%
        h.scale             = fread(fid,1,'double');
        h.veogcorrect       = fread(fid,1,'uint8');
        h.heogcorrect       = fread(fid,1,'uint8');
        h.aux1correct       = fread(fid,1,'uint8');
        h.aux2correct       = fread(fid,1,'uint8');
        h.veogtrig          = fread(fid,1,'float');
        h.heogtrig          = fread(fid,1,'float');
        h.aux1trig          = fread(fid,1,'float');
        h.aux2trig          = fread(fid,1,'float');
        h.heogchnl          = fread(fid,1,'short');
        h.veogchnl          = fread(fid,1,'short');
        h.aux1chnl          = fread(fid,1,'short');
        h.aux2chnl          = fread(fid,1,'short');
        h.veogdir           = fread(fid,1,'uint8');
        h.heogdir           = fread(fid,1,'uint8');
        h.aux1dir           = fread(fid,1,'uint8');
        h.aux2dir           = fread(fid,1,'uint8');
        h.veog_n            = fread(fid,1,'short');
        h.heog_n            = fread(fid,1,'short');
        h.aux1_n            = fread(fid,1,'short');
        h.aux2_n            = fread(fid,1,'short');
        h.veogmaxcnt        = fread(fid,1,'short');
        h.heogmaxcnt        = fread(fid,1,'short');
        h.aux1maxcnt        = fread(fid,1,'short');
        h.aux2maxcnt        = fread(fid,1,'short');
        h.veogmethod        = fread(fid,1,'uint8');
        h.heogmethod        = fread(fid,1,'uint8');
        h.aux1method        = fread(fid,1,'uint8');
        h.aux2method        = fread(fid,1,'uint8');
        h.ampsensitivity    = fread(fid,1,'float');
        h.lowpass           = fread(fid,1,'uint8');	%%%
        h.highpass          = fread(fid,1,'uint8');	%%%
        h.notch             = fread(fid,1,'uint8');	%%%
        h.autoclipadd       = fread(fid,1,'uint8');
        h.baseline          = fread(fid,1,'uint8');	%%%
        h.offstart          = fread(fid,1,'float');
        h.offstop           = fread(fid,1,'float');
        h.reject            = fread(fid,1,'uint8');
        h.rejstart          = fread(fid,1,'float');
        h.rejstop           = fread(fid,1,'float');
        h.rejmin            = fread(fid,1,'float');
        h.rejmax            = fread(fid,1,'float');
        h.trigtype          = fread(fid,1,'uint8');
        h.trigval           = fread(fid,1,'float');
        h.trigchnl          = fread(fid,1,'uint8');	%%%
        h.trigmask          = fread(fid,1,'short');
        h.trigisi           = fread(fid,1,'float');
        h.trigmin           = fread(fid,1,'float');
        h.trigmax           = fread(fid,1,'float');
        h.trigdir           = fread(fid,1,'uint8');
        h.autoscale         = fread(fid,1,'uint8');
        h.n2                = fread(fid,1,'short');
        h.dir               = fread(fid,1,'uint8');
        h.dispmin           = fread(fid,1,'float');
        h.dispmax           = fread(fid,1,'float');
        h.xmin              = fread(fid,1,'float');
        h.xmax              = fread(fid,1,'float');
        h.automin           = fread(fid,1,'float');
        h.automax           = fread(fid,1,'float');
        h.zmin              = fread(fid,1,'float');
        h.zmax              = fread(fid,1,'float');
        h.lowcut            = fread(fid,1,'float');	%%%
        h.highcut           = fread(fid,1,'float');	%%%
        h.common            = fread(fid,1,'uint8');
        h.savemode          = fread(fid,1,'uint8');
        h.manmode           = fread(fid,1,'uint8');
        h.ref               = fread(fid,10,'uint8');
        h.rectify           = fread(fid,1,'uint8');
        h.displayxmin       = fread(fid,1,'float');
        h.displayxmax       = fread(fid,1,'float');
        h.phase             = fread(fid,1,'uint8');
        h.screen            = fread(fid,16,'uint8');
        h.calmode           = fread(fid,1,'short');
        h.calmethod         = fread(fid,1,'short');
        h.calupdate         = fread(fid,1,'short');
        h.calbaseline       = fread(fid,1,'short');
        h.calsweeps         = fread(fid,1,'short');
        h.calattenuator     = fread(fid,1,'float');
        h.calpulsevolt      = fread(fid,1,'float');
        h.calpulsestart     = fread(fid,1,'float');
        h.calpulsestop      = fread(fid,1,'float');
        h.calfreq           = fread(fid,1,'float');
        h.taskfile          = fread(fid,34,'uint8');
        h.seqfile           = fread(fid,34,'uint8');
        h.spectmethod       = fread(fid,1,'uint8');
        h.spectscaling      = fread(fid,1,'uint8');
        h.spectwindow       = fread(fid,1,'uint8');
        h.spectwinlength    = fread(fid,1,'float');
        h.spectorder        = fread(fid,1,'uint8');
        h.notchfilter       = fread(fid,1,'uint8');	%%%
        h.headgain          = fread(fid,1,'short');	%%%	
        h.additionalfiles   = fread(fid,1,'int');
        h.unused            = fread(fid,5,'uint8');
        h.fspstopmethod     = fread(fid,1,'short');
        h.fspstopmode       = fread(fid,1,'short');
        h.fspfvalue         = fread(fid,1,'float');
        h.fsppoint          = fread(fid,1,'short');
        h.fspblocksize      = fread(fid,1,'short');
        h.fspp1             = fread(fid,1,'ushort');
        h.fspp2             = fread(fid,1,'ushort');
        h.fspalpha          = fread(fid,1,'float');
        h.fspnoise          = fread(fid,1,'float');
        h.fspv1             = fread(fid,1,'short');
        h.montage           = fread(fid,40,'uint8');
        h.eventfile         = fread(fid,40,'uint8');
        h.fratio            = fread(fid,1,'float');
        h.minor_rev         = fread(fid,1,'uint8');	%%%
        h.eegupdate         = fread(fid,1,'short');
        h.compressed        = fread(fid,1,'uint8');
        h.xscale            = fread(fid,1,'float');
        h.yscale            = fread(fid,1,'float');
        h.xsize             = fread(fid,1,'float');
        h.ysize             = fread(fid,1,'float');
        h.acmode            = fread(fid,1,'uint8');
        h.commonchnl        = fread(fid,1,'uint8');
        h.xtics             = fread(fid,1,'uint8');
        h.xrange            = fread(fid,1,'uint8');
        h.ytics             = fread(fid,1,'uint8');
        h.yrange            = fread(fid,1,'uint8');
        h.xscalevalue       = fread(fid,1,'float');
        h.xscaleinterval    = fread(fid,1,'float');
        h.yscalevalue       = fread(fid,1,'float');
        h.yscaleinterval    = fread(fid,1,'float');
        h.scaletoolx1       = fread(fid,1,'float');
        h.scaletooly1       = fread(fid,1,'float');
        h.scaletoolx2       = fread(fid,1,'float');
        h.scaletooly2       = fread(fid,1,'float');
        h.port              = fread(fid,1,'short');
%        h.numsamples        = fread(fid,1,'uint32');	%%%
        h.numsamples        = fread(fid,1,'float32');	%%%

        h.filterflag        = fread(fid,1,'uint8');	%%%
        h.lowcutoff         = fread(fid,1,'float');	%%%
        h.lowpoles          = fread(fid,1,'short');	
        h.highcutoff        = fread(fid,1,'float');	%%%
        h.highpoles         = fread(fid,1,'short');
        h.filtertype        = fread(fid,1,'uint8');
        h.filterdomain      = fread(fid,1,'uint8');
        h.snrflag           = fread(fid,1,'uint8');
        h.coherenceflag     = fread(fid,1,'uint8');
        h.continuoustype    = fread(fid,1,'uint8');
        h.eventtablepos     = fread(fid,1,'int32');	%%%
        h.continuousseconds = fread(fid,1,'float');
        h.channeloffset     = fread(fid,1,'uint32');
        h.autocorrectflag   = fread(fid,1,'uint8');
        h.dcthreshold       = fread(fid,1,'uint8');
        
        if ftell(fid)~=900,
                warning(['supicous Neuroscan file ',FILENAME]);
        end;
        
        for n = 1:h.nchannels,%h.nchannels
                e.lab(1:10,n)         = fread(fid,10,'uint8');
                e.reference(1,n)      = fread(fid,1,'uint8');
                e.skip(1,n)           = fread(fid,1,'uint8');
                e.reject(1,n)         = fread(fid,1,'uint8');
                e.display(1,n)        = fread(fid,1,'uint8');
                e.bad(1,n)            = fread(fid,1,'uint8');
                e.n(1,n)              = fread(fid,1,'ushort');
                e.avg_reference(1,n)  = fread(fid,1,'uint8');
                e.clipadd(1,n)        = fread(fid,1,'uint8');
                e.x_coord(1,n)        = fread(fid,1,'float');
                e.y_coord(1,n)        = fread(fid,1,'float');
                e.veog_wt(1,n)        = fread(fid,1,'float');
                e.veog_std(1,n)       = fread(fid,1,'float');
                e.snr(1,n)            = fread(fid,1,'float');
                e.heog_wt(1,n)        = fread(fid,1,'float');
                e.heog_std(1,n)       = fread(fid,1,'float');
                e.baseline(1,n)       = fread(fid,1,'short');
                e.filtered(1,n)       = fread(fid,1,'uint8');
                e.fsp(1,n)            = fread(fid,1,'uint8');
                e.aux1_wt(1,n)        = fread(fid,1,'float');
                e.aux1_std(1,n)       = fread(fid,1,'float');
                e.sensitivity(1,n)    = fread(fid,1,'float');
                e.gain(1,n)           = fread(fid,1,'uint8');
                e.hipass(1,n)         = fread(fid,1,'uint8');
                e.lopass(1,n)         = fread(fid,1,'uint8');
                e.page(1,n)           = fread(fid,1,'uint8');
                e.size(1,n)           = fread(fid,1,'uint8');
                e.impedance(1,n)      = fread(fid,1,'uint8');
                e.physicalchnl(1,n)   = fread(fid,1,'uint8');
                e.rectify(1,n)        = fread(fid,1,'uint8');
                e.calib(1,n)          = fread(fid,1,'float');
        end
        
        if ftell(fid)~=(900+h.nchannels*75),	
	% this check does not work in the currenct CVS-version of Octave	
                warning(['supicous Neuroscan file ',FILENAME]);
        end;
        
end;

CNT.VERSION = str2double(char(h.rev(8:12)'));
CNT.CNT.type = h.type;
CNT.PID = h.id;
CNT.ID.Operator = char(h.oper');	%
CNT.ID.Doctor   = char(h.doctor');	%
CNT.ID.referral = char(h.referral');	%
CNT.ID.Hospital = char(h.hospital');	%
CNT.Patient.Name= char(h.patient');	%
CNT.Patient.Age = h.age;	%
CNT.Patient.Sex = char(h.sex');	%
CNT.Patient.Handedness=char(h.hand');	%	
CNT.Patient.Medication=char(h.med');	%	
CNT.Patient.Classification=char(h.category');	%	 
CNT.Patient.State=char(h.state');	%	
CNT.Session.Label=char(h.label');	%	
CNT.Date=char(h.date');	%	
CNT.Time=char(h.time');	%	
if any(CNT.Date=='/') & ~any(CNT.Date=='-') & ~any(CNT.Date=='.'),	%% data format MM/DD/YYYY
	CNT.T0 = [str2double(CNT.Date(7:length(CNT.Date))),str2double(CNT.Date(1:2)),str2double(CNT.Date(4:5)),str2double(CNT.Time(1:2)),str2double(CNT.Time(4:5)),str2double(CNT.Time(7:8))];
else %if ~any(CNT.Date=='/'),	%% data format DD-MM-YYYY or DD.MM.YYYY 
	CNT.T0 = [str2double(CNT.Date(7:length(CNT.Date))),str2double(CNT.Date(4:5)),str2double(CNT.Date(1:2)),str2double(CNT.Time(1:2)),str2double(CNT.Time(4:5)),str2double(CNT.Time(7:8))];
end;
% check year
if     CNT.T0(1) > 99,
elseif CNT.T0(1) > 80, 	CNT.T0(1) = CNT.T0(1) + 1900;
else			CNT.T0(1) = CNT.T0(1) + 2000;
end;
% check day & month
if CNT.T0(2)>12,
        fprintf(2, 'Warning CNTOPEN: month and day were mixed up %i-%i-%i-%i-%i-%i \n',CNT.T0);
        CNT.T0(2:3) = CNT.T0([3,2]);
end;

CNT.NS = h.nchannels;	%	
CNT.SampleRate=h.rate;	% D-to-A rate	
CNT.Scale=h.scale;	% scale factor for calibration
CNT.Scale2=h.ampsensitivity;
CNT.HeadLen = 900 + 75*CNT.NS;
%CNT.PhysDim = repmat({'µV'},CNT.NS,1);
CNT.PhysDimCode = repmat(4275,CNT.NS,1); %% uV 

% Scan4.3->Edit->Overall Setup->Amplifier->Notch->Off/50Hz/60Hz
tmp = [0,50,60]; 
CNT.Filter.Notch  = tmp(h.notchfilter+1);

% Scan4.3->Edit->Overall Setup->Amplifier->AC/DC
CNT.Filter.ACmode = h.acmode;

% Scan4.3->Edit->Overall Setup->Amplifier->DC Auto Correction
CNT.Filter.DCauto = h.autocorrectflag;

% Scan4.3->Edit->Overall Setup->Amplifier->Amplifier Settings->Low Pass
tmp = [30, 40, 50, 70, 100, 200, 500, 1000, 1500, 2000, 2500, 3000]; % LOWPASS
CNT.Filter.LowPass = tmp(e.lopass+1);

CNT.CNT.Filter.LowPass = tmp(h.lowpass+1); % ???

% Scan4.3->Edit->Overall Setup->Amplifier->Amplifier Settings->High Pass
tmp = [NaN, 0, .05, .1, .15, .3, 1, 5, 10, 30, 100, 150, 300]; %HIGHPASS
CNT.Filter.HighPass = tmp(e.hipass+1);
if h.acmode,
        CNT.Filter.HighPass(e.hipass==0) = .05;
end;

CNT.CNT.Filter.HighPass = tmp(h.highpass+1); % ???

 % ???
CNT.CNT.Filter.LowCutOff  = h.lowcutoff;
CNT.CNT.Filter.HighCutOff = h.highcutoff;
CNT.CNT.Filter.NotchOn = h.filterflag;
CNT.CNT.Filter.ON   = [e(:).filtered];
CNT.CNT.minor_revision = h.minor_rev;
CNT.CNT.EventTablePos  = h.eventtablepos;

CNT.Label = char(e.lab');

CNT.FILE.POS = 0;
if strcmp(upper(CNT.FILE.Ext),'AVG'),
        if (h.type~=0),
		fprintf(2,'Warning CNTOPEN: filetype %i does not match file extension (%s).\n',h.type,CNT.FILE.Ext); 
	end;
	CNT.TYPE='AVG';
        CNT.AS.endpos = 1;
	CNT.NRec = 1;
        CNT.SPR  = h.pnts;
        CNT.Cal  = e.calib./e.n;   % scaling
        CNT.Calib = sparse(2:CNT.NS+1,1:CNT.NS,CNT.Cal);
	CNT.AS.bpb = h.pnts*h.nchannels*4+5;
	CNT.AS.spb = h.pnts*h.nchannels;
	CNT.Dur = CNT.SPR/CNT.SampleRate;
        CNT.GDFTYP = 16; %'float32';
        
elseif strcmp(upper(CNT.FILE.Ext),'COH')        
        warning('.COH data not supported yet')
        CNT.COH.directory = fread(CNT.FILE.FID,[CNT.NS,CNT.NS],'int32');
        CNT.SPR = h.pnts;
        CNT.GDFTYP = 16; %'float32';
        
elseif strcmp(upper(CNT.FILE.Ext),'CSA')        
        warning('.CSA data not supported yet')
        CNT.SPR  = h.pnts;
        CNT.NRec = h.compsweeps;
        CNT.GDFTYP = 16; %'float32';
        
elseif strcmp(upper(CNT.FILE.Ext),'EEG'),
	if (h.type~=1),
		fprintf(2,'Warning CNTOPEN: filetype %i does not match file extension (%s).\n',h.type,CNT.FILE.Ext); 
	end;
	CNT.TYPE   = 'EEG';
        CNT.SPR    = h.pnts;
        CNT.NRec   = h.compsweeps;
        CNT.AS.spb = CNT.NS*CNT.SPR;	% Samples per Block
        
        % Sometimes h.eventtablepos seems to need a correction, also I've not figured out why. 
        % The Manual SCAN 4.2 Vol II, Page Headers-7 refers to "286 SCAN manual". Maybe this could bring a clarification. 
        % Anyway, the following code deals with the problem.   
        CNT.AS.bpb = -1;

        if CNT.CNT.minor_revision==12,
                CNT.AS.bpb = 2*CNT.AS.spb+1+2+2+4+2+2;
                CNT.GDFTYP = 3; %'int16';
                % correct(?) eventtablepos
                h.eventtablepos = CNT.HeadLen + CNT.NRec*CNT.AS.bpb;    	    
        else
                if CNT.CNT.minor_revision~=16,
                        fprintf(CNT.FILE.stderr,'Warning CNTOPEN: EEG-Format Minor-Revision %i not tested.\n',CNT.CNT.minor_revision);
                end;
		
                tmp = (CNT.AS.spb*2+1+2+2+4+2+2);
                if (h.eventtablepos-CNT.HeadLen)==(tmp*CNT.NRec),
                        CNT.AS.bpb = tmp;
                        CNT.GDFTYP = 3; %'int16';
                end;
		
                tmp = (CNT.AS.spb*4+1+2+2+4+2+2);
                if (h.eventtablepos-CNT.HeadLen)==(tmp*CNT.NRec),
                        CNT.AS.bpb = tmp;
                        CNT.GDFTYP = 5; %'int32';
                end;
        end; 
        if CNT.AS.bpb < 0;
                fprintf(CNT.FILE.stderr,'Error CNTOPEN: header information of file %s corrupted.\n',CNT.FileName);
                fclose(CNT.FILE.FID);
                CNT.FILE.FID = -1;
                return;
	end;
        
        CNT.Calib = [-[e.baseline];eye(CNT.NS)]*diag([e.sensitivity].*[e.calib]/204.8);
        CNT.AS.endpos = CNT.NRec;
	CNT.FLAG.TRIGGERED = 1;
	CNT.Dur = CNT.SPR/CNT.SampleRate;
        
elseif  strcmp(upper(CNT.FILE.Ext),'CNT'),
        CNT.TYPE = 'CNT';
        CNT.SPR   = h.numsamples;
        %CNT.SPR    = h.pnts;
        CNT.NRec   = h.compsweeps;
        %disp([h.eventtablepos,CNT.HeadLen,CNT.NS,h.pnts,CNT.NRec,CNT.SampleRate,h.type,CNT.CNT.minor_revision])
        %disp([CNT.NS,h.pnts,h.compsweeps,h.numsamples,h.type,CNT.CNT.minor_revision])
        CNT.CNT.h = h; 
        
        if (CNT.CNT.minor_revision==8),
	        CNT.GDFTYP = 3; %'int16';
                h.numsamples; % might have some meaning 
        elseif (CNT.CNT.minor_revision==12),
                CNT.GDFTYP = 3; %'int16';
        elseif (CNT.CNT.minor_revision==16),
                CNT.GDFTYP = 3; %'int16';
        else
                fprintf(CNT.FILE.stderr,'Warning CNTOPEN: EEG-Format Minor-Revision %i not tested.\n',CNT.CNT.minor_revision);
        end;

        if ~isempty(strfind(arg2,'32')),        % force 32 bit
                % if h.type==184
                CNT.GDFTYP = 5; %'int32';
                CNT.AS.bpb = CNT.NS*4;	% Bytes per Block
        else
                CNT.GDFTYP = 3; %'int16';
                CNT.AS.bpb = CNT.NS*2;	% Bytes per Block
        end;

        if h.eventtablepos>CNT.FILE.size,
                warning('CNTOPEN: eventtablepos %i after end of file %i: \n',h.eventtablepos,CNT.FILE.size);
        end; 
        
	CNT.AS.spb = CNT.NS;	% Samples per Block
	CNT.AS.EVENTTABLEPOS = h.eventtablepos;
	if h.eventtablepos>CNT.FILE.size,
                fprintf(CNT.FILE.stderr,'Warning CNTOPEN: %s is corrupted:\n    - position of eventtable (%i) past end of file (%i).\n',CNT.FileName,h.eventtablepos,CNT.FILE.size);
	        CNT.SPR    = floor((CNT.FILE.size-CNT.HeadLen)/CNT.AS.bpb);
        else
	        CNT.SPR    = (h.eventtablepos-CNT.HeadLen)/CNT.AS.bpb;
	end;	
	CNT.AS.endpos = CNT.SPR;
        
        CNT.NRec   = 1;
	CNT.Calib  = [-[e.baseline];eye(CNT.NS)]*diag([e.sensitivity].*[e.calib]/204.8);
	CNT.FLAG.TRIGGERED = 0;	        
	CNT.Dur    = 1/CNT.SampleRate;

	if all(CNT.GDFTYP==3)
		CNT.DigMax = repmat(32767,1,CNT.NS); 
		CNT.DigMin = repmat(-32768,1,CNT.NS); 
	else
		CNT.DigMax = repmat(2^23-1,1,CNT.NS); 
		CNT.DigMin = repmat(-2^23,1,CNT.NS); 
	end; 
	CNT.PhysMax = [1,CNT.DigMax]*CNT.Calib; 	
	CNT.PhysMin = [1,CNT.DigMin]*CNT.Calib; 	

        %%%%% read event table 
        CNT.EVENT.TYP = [];
        CNT.EVENT.POS = [];

	CNT.EVENT.TeegSize = 0;
        status = fseek(CNT.FILE.FID,h.eventtablepos,'bof');
        if ~status,
                [CNT.EVENT.TeegType,c1] = fread(fid,1,'uint8');		
                [CNT.EVENT.TeegSize,c2] = fread(fid,1,'int32');	
                [CNT.EVENT.TeegOffset,c3] = fread(fid,1,'int32');
	end;	
	
        k = 0; 
        K = 1;
        TEEG = [];
        while (K < CNT.EVENT.TeegSize),
		k = k+1;
                Teeg.Stimtype = fread(fid,1,'int16');        
                Teeg.Keyboard = fread(fid,1,'uint8');        
                tmp = fread(fid,1,'uint8');
                Teeg.tmp = tmp;         
                Teeg.KeyPad = rem(tmp,16); %bitand(tmp,15);
                Teeg.Accept = (fix(tmp/16))==13; % (bitshift(tmp,-4)==13);  % 0xd = accept, 0xc = reject 
                        
                Teeg.Offset = fread(fid,1,'int32');        
                K = K + 8;
                if any(CNT.EVENT.TeegType==[2:3]),
                        Teeg.Type       =  fread(fid,1,'int16');        
                        Teeg.Code       =  fread(fid,1,'int16');        
                        Teeg.Latency    =  fread(fid,1,'float32');        
                        Teeg.EpochEvent =  fread(fid,1,'uint8');        
                        Teeg.Accept2    =  fread(fid,1,'uint8');        
                        Teeg.Accuracy   =  fread(fid,1,'uint8');        
                        K = K + 11;
                end;    
		if k==1,
	                TEEG = Teeg;
    		else
		        TEEG(k) = Teeg;
		end;
        end;
	
        if length(TEEG) > 0,
                CNT.EVENT.TEEG = TEEG'; 
                %CNT.EVENT.TYP = [TEEG(:).Stimtype]';
                stim = [TEEG(:).Stimtype]';
                resp = [TEEG(:).KeyPad]';
                ix   = find((stim>0) & (resp>0));
                if ~isempty(ix)
                	fprintf(CNT.FILE.stderr,'Warning SOPEN(CNT): in some events, both response and stimulus are non zero. Response code ignored.\n');
                	ix',
                end;	
                
                CNT.EVENT.TYP  = stim + (resp+128*(stim==0));
                if CNT.EVENT.TeegType==3,
                        CNT.EVENT.POS =  [TEEG(:).Offset]';
                else
                        CNT.EVENT.POS = ([TEEG(:).Offset]' - CNT.HeadLen) ./ CNT.AS.bpb;
                end;
        end;
end;

% set file pointer to the beginning of the data block
status = fseek(CNT.FILE.FID, CNT.HeadLen, 'bof');
if status,
        fprintf(CNT.FILE.stderr,'Warning CNTOPEN: I/O error in file %s\n',CNT.FileName);
end;        
CNT.h = h;