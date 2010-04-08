function cdf = spm_read_netcdf(fname)
% Read the header information from a NetCDF file into a data structure.
% FORMAT cdf = spm_read_netcdf(fname)
% fname - name of NetCDF file
% cdf   - data structure
%
% See: http://www.unidata.ucar.edu/packages/netcdf/
% _______________________________________________________________________
% @(#)spm_read_netcdf.m	2.1 John Ashburner 02/07/30

dsiz     = [1 1 2 4 4 8];
fp=fopen(fname,'r','ieee-be');
if fp==-1,
	cdf = [];
	return;
end;

% Return null if not a CDF file.
%-----------------------------------------------------------------------
mgc = fread(fp,4,'uchar')';
if ~all(['CDF' 1] == mgc),
	cdf = [];
	fclose(fp);
	return;
end

% I've no idea what this is for
numrecs = fread(fp,1,'uint32');

cdf = struct('numrecs',numrecs,'dim_array',[], 'gatt_array',[], 'var_array', []);

dt = fread(fp,1,'uint32');
if dt == 10,
	% Dimensions
	nelem = fread(fp,1,'uint32');
	for j=1:nelem,
		str   = readname(fp);
		dim_length = fread(fp,1,'uint32');
		cdf.dim_array(j).name       = str;
		cdf.dim_array(j).dim_length = dim_length;
	end;
	dt   = fread(fp,1,'uint32');
end

while ~dt, dt   = fread(fp,1,'uint32'); end;

if dt == 12,
	% Attributes
	nelem = fread(fp,1,'uint32');
	for j=1:nelem,
		str    = readname(fp);
		nc_type= fread(fp,1,'uint32');
		nnelem = fread(fp,1,'uint32');
		val    = fread(fp,nnelem,dtypestr(nc_type));
		if nc_type == 2, val = deblank([val' ' ']); end
		padding= fread(fp,ceil(nnelem*dsiz(nc_type)/4)*4-nnelem*dsiz(nc_type),'uchar');
		cdf.gatt_array(j).name    = str;
		cdf.gatt_array(j).nc_type = nc_type;
		cdf.gatt_array(j).val     = val;
	end;
	dt   = fread(fp,1,'uint32');
end

while ~dt, dt   = fread(fp,1,'uint32'); end;

if dt == 11,
	% Variables
	nelem = fread(fp,1,'uint32');
	for j=1:nelem,
		str    = readname(fp);
		nnelem = fread(fp,1,'uint32');
		val    = fread(fp,nnelem,'uint32');
		cdf.var_array(j).name    = str;
		cdf.var_array(j).dimid   = val+1;
		cdf.var_array(j).nc_type = 0;
		cdf.var_array(j).vsize   = 0;
		cdf.var_array(j).begin   = 0;
		dt0    = fread(fp,1,'uint32');
		if dt0 == 12,
			nelem0 = fread(fp,1,'uint32');
			for jj=1:nelem0,
				str    = readname(fp);
				nc_type= fread(fp,1,'uint32');
				nnelem = fread(fp,1,'uint32');
				val    = fread(fp,nnelem,dtypestr(nc_type));
				if nc_type == 2, val = deblank([val' ' ']); end
				padding= fread(fp,...
					ceil(nnelem*dsiz(nc_type)/4)*4-nnelem*dsiz(nc_type),'uchar');
				cdf.var_array(j).vatt_array(jj).name    = str;
				cdf.var_array(j).vatt_array(jj).nc_type = nc_type;
				cdf.var_array(j).vatt_array(jj).val     = val;
			end;
			dt0    = fread(fp,1,'uint32');
		end;
		cdf.var_array(j).nc_type  = dt0;
		cdf.var_array(j).vsize = fread(fp,1,'uint32');
		cdf.var_array(j).begin = fread(fp,1,'uint32');
	end;
	dt   = fread(fp,1,'uint32');
end;

fclose(fp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function str = dtypestr(i)
% Returns a string appropriate for reading or writing the CDF data-type.
types = str2mat('uint8','uint8','int16','int32','float','double');
str   = deblank(types(i,:));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function name = readname(fp)
% Extracts a name from a CDF file pointed to at the right location by
% fp.
stlen  = fread(fp,1,'uint32');
name   = deblank([fread(fp,stlen,'uchar')' ' ']);
padding= fread(fp,ceil(stlen/4)*4-stlen,'uchar');
return;
%_______________________________________________________________________
