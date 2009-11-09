function [hdr,otherendian] = spm_read_hdr(fname)
% Read (SPM customised) Analyze header
% FORMAT [hdr,otherendian] = spm_read_hdr(fname)
% fname       - .hdr filename
% hdr         - structure containing Analyze header
% otherendian - byte swapping necessary flag
%_______________________________________________________________________
% @(#)spm_read_hdr.m	2.2 John Ashburner 03/07/17

fid         = fopen(fname,'r','native');
otherendian = 0;
if (fid > 0)
	dime = read_dime(fid);
	if dime.dim(1)<0 | dime.dim(1)>15, % Appears to be other-endian
		% Re-open other-endian
		fclose(fid);
		if spm_platform('bigend'), fid = fopen(fname,'r','ieee-le');
		else,                      fid = fopen(fname,'r','ieee-be'); end;
		otherendian = 1;
		dime = read_dime(fid);
	end;
	hk       = read_hk(fid);
	hist     = read_hist(fid);
	hdr.hk   = hk;
	hdr.dime = dime;
	hdr.hist = hist;

	% SPM specific bit - unused
	%if hdr.hk.sizeof_hdr > 348,
	%	spmf = read_spmf(fid,dime.dim(5));
	%	if ~isempty(spmf),
	%		hdr.spmf = spmf;
	%	end;
	%end;

	fclose(fid);
else,
	hdr = [];
	otherendian = NaN;
	%error(['Problem opening header file (' fopen(fid) ').']);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function hk = read_hk(fid)
% read (struct) header_key
%-----------------------------------------------------------------------
fseek(fid,0,'bof');
hk.sizeof_hdr 		= fread(fid,1,'int32');
hk.data_type  		= mysetstr(fread(fid,10,'uchar'))';
hk.db_name    		= mysetstr(fread(fid,18,'uchar'))';
hk.extents    		= fread(fid,1,'int32');
hk.session_error	= fread(fid,1,'int16');
hk.regular			= mysetstr(fread(fid,1,'uchar'))';
hk.hkey_un0			= mysetstr(fread(fid,1,'uchar'))';
if isempty(hk.hkey_un0), error(['Problem reading "hk" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function dime = read_dime(fid)
% read (struct) image_dimension
%-----------------------------------------------------------------------
fseek(fid,40,'bof');
dime.dim		= fread(fid,8,'int16')';
dime.vox_units	= mysetstr(fread(fid,4,'uchar'))';
dime.cal_units	= mysetstr(fread(fid,8,'uchar'))';
dime.unused1	= fread(fid,1,'int16');
dime.datatype	= fread(fid,1,'int16');
dime.bitpix		= fread(fid,1,'int16');
dime.dim_un0	= fread(fid,1,'int16');
dime.pixdim		= fread(fid,8,'float')';
dime.vox_offset	= fread(fid,1,'float');
dime.funused1	= fread(fid,1,'float');
dime.funused2	= fread(fid,1,'float');
dime.funused3	= fread(fid,1,'float');
dime.cal_max	= fread(fid,1,'float');
dime.cal_min	= fread(fid,1,'float');
dime.compressed	= fread(fid,1,'int32');
dime.verified	= fread(fid,1,'int32');
dime.glmax		= fread(fid,1,'int32');
dime.glmin		= fread(fid,1,'int32');
if isempty(dime.glmin), error(['Problem reading "dime" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function hist = read_hist(fid)
% read (struct) data_history
%-----------------------------------------------------------------------
fseek(fid,148,'bof');
hist.descrip	= mysetstr(fread(fid,80,'uchar'))';
hist.aux_file	= mysetstr(fread(fid,24,'uchar'))';
hist.orient		= fread(fid,1,'uchar');
hist.origin		= fread(fid,5,'int16')';
hist.generated	= mysetstr(fread(fid,10,'uchar'))';
hist.scannum	= mysetstr(fread(fid,10,'uchar'))';
hist.patient_id	= mysetstr(fread(fid,10,'uchar'))';
hist.exp_date	= mysetstr(fread(fid,10,'uchar'))';
hist.exp_time	= mysetstr(fread(fid,10,'uchar'))';
hist.hist_un0	= mysetstr(fread(fid,3,'uchar'))';
hist.views		= fread(fid,1,'int32');
hist.vols_added	= fread(fid,1,'int32');
hist.start_field= fread(fid,1,'int32');
hist.field_skip	= fread(fid,1,'int32');
hist.omax		= fread(fid,1,'int32');
hist.omin		= fread(fid,1,'int32');
hist.smax		= fread(fid,1,'int32');
hist.smin		= fread(fid,1,'int32');
if isempty(hist.smin), error(['Problem reading "hist" of header file (' fopen(fid) ').']); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function spmf = read_spmf(fid,n)
% Read SPM specific fields
% This bit may be used in the future for extending the Analyze header.

fseek(fid,348,'bof');
mgc = fread(fid,1,'int32');    % Magic number
if mgc ~= 20020417, spmf = []; return; end;

for j=1:n,
	spmf(j).mat    = fread(fid,16,'double'); % Orientation information
	spmf(j).unused = fread(fid,384,'uchar'); % Extra unused stuff
	if length(spmf(j).unused)<384,
		error(['Problem reading "spmf" of header file (' fopen(fid) ').']);
	end;
 	spmf(j).mat = reshape(spmf(j).mat,[4 4]);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function out = mysetstr(in)
tmp = find(in == 0);
tmp = min([min(tmp) length(in)]);
out = setstr([in(1:tmp)' zeros(1,length(in)-(tmp))])';
return;
%_______________________________________________________________________
%_______________________________________________________________________
