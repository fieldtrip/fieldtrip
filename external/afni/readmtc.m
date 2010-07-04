function mtc = readmtc (filename, o1d)
% mtc = readmtc (filename, o1d);
% Read BrainVoyager's mesh time course (mtc) files
% filename: name of .mtc file
% o1d: (optional), name of output .1D (AFNI's text format output) version
% mtc: output stucture with contents of filename
%
% By ZSS, SSCC/NIMH/NIH 
% based on Hester's positing at
%
% http://www.brainvoyager.com/ubb/Forum8/HTML/000059.html
%
%The format of MTC files
%
%Mesh Time Course (*.mtc) files are generated from BrainVoyager QX version 
%1.0. The files are saved in little endian byte order.
%
%MTC header
%
%BYTES : DATA TYPE : (DEFAULT) : DESCRIPTION
%
%4 : int: (1) : version number
%4 : int: number of vertices
%4 : int: number of volumes
%x * 1 : char : name of source VTC file' (ends with '0')
%x * 1 : char : name of linked protocol file; if not available: ''; (ends with '0')
%4 : int : hemodynamic delay
%4 : float : TR - Repetition Time
%4 : float : delta parameter for hemodynamic response function
%4 : float : tau parameter for hemodynamic response function
%4 : int : segment size (intervals per stimulus block/event)
%4 : int : segment offset - first datapoint used with 
%protocol
%1 : char : (1) : datatype of MTC data; 1 = 'float'
%
%MTC data
%
%The data are organized per vertex: all datapoints per vertex are saved 
%together, so time runs fastest (has the smallest loop).
%
%Nr of vertices
%
%Nr of timepoints
%
%***
%
%The *.mtc file format description is also added to the updated file format documentation as available on our ftp server.
%
%[This message has been edited by Hester Breman (edited 29 October 2004).]
%

if (nargin == 1),
   o1d = ''
end

if (~isempty(o1d)),
   o1d = Remove1DExtension(o1d);
   o1d = sprintf('%s.1D.dset', o1d);
   if (o1d(1) ~= filesep), o1d = sprintf('.%c%s', filesep, o1d); end
   if (filexist(o1d)),
      fprintf(2,'Error: Output file %s exists.\nWill not overwrite.\n', o1d);
      return;
   end
end

fprintf(2,'Reading from %s...\n', filename);

to = clock;
fp = fopen(filename,'r');
if (fp == -1) 
	fprintf(1,'\nError opening %s\n',filename);
	return;
end


%% read some header fields 
mtc.version = fread(fp,1,'int32',0,'ieee-le');
mtc.numvert = fread(fp,1,'int32',0,'ieee-le');
mtc.numvol = fread(fp,1,'int32',0,'ieee-le');
cnt = 1;
mtc.vtcname = '';
c=fread(fp,1,'char',0,'ieee-le');
%fprintf(2,'%d-->%c\n',c,c);
while ( c ~= 0 ),
   mtc.vtcname(cnt) = c;
   c=fread(fp,1,'char',0,'ieee-le');
   %fprintf(2,'%d-->%c\n',c,c);
   cnt = cnt + 1;
end
cnt = 1;
mtc.lpfname = '';
c=fread(fp,1,'char',0,'ieee-le');
%fprintf(2,'%d-->%c\n',c,c);
while ( c ~= 0 ),
   mtc.lpfname(cnt) = c;
   c=fread(fp,1,'char',0,'ieee-le');
   %fprintf(2,'%d-->%c\n',c,c);
  cnt = cnt + 1;
end

mtc.hemdelay = fread(fp,1,'int32',0,'ieee-le');
mtc.tr = fread(fp,1,'float32',0,'ieee-le');
mtc.hrfdelta = fread(fp,1,'float32',0,'ieee-le');
mtc.hrftau = fread(fp,1,'float32',0,'ieee-le');
mtc.segsize = fread(fp,1,'int32',0,'ieee-le');
mtc.sefoffset = fread(fp,1,'int32',0,'ieee-le'); 
mtc.datatype = fread(fp,1,'char',0,'ieee-le');

if (mtc.datatype ~= 1),
   fprintf(2,'Datatype unexpected.\n');
   return;
end
mtc.data = zeros(mtc.numvert, mtc.numvol);
for (i=1:1:mtc.numvert),
   mtc.data(i,:) = fread(fp,mtc.numvol,'float32',0,'ieee-le');
end
fclose(fp);

[h,m,s] = sectotime(etime(clock, to));
fprintf(2,'Reading took %g:%g:%g.\n', h,m,s);

if (~isempty(o1d)),
   fprintf(2,'Writing to %s ...\n', o1d);
   to = clock;
   fo = fopen(o1d,'w');
   fprintf(fo,'#BrainVoyager .mtc file transformed to .1D format with readmtc.m\n');
   fprintf(fo,'#version: %g\n', mtc.version);
   fprintf(fo,'#numvert: %g\n', mtc.numvert);
   fprintf(fo,'#numvol: %g\n', mtc.numvol);
   fprintf(fo,'#vtcname: %s\n', mtc.vtcname);
   fprintf(fo,'#lpfname: %s\n', mtc.lpfname);
   fprintf(fo,'#hemdelay: %g\n', mtc.hemdelay);
   fprintf(fo,'#tr: %g\n', mtc.tr);
   fprintf(fo,'#hrfdelta: %g\n', mtc.hrfdelta);
   fprintf(fo,'#hrftau: %g\n', mtc.hrftau);
   fprintf(fo,'#segsize: %g\n', mtc.segsize);
   fprintf(fo,'#sefoffset: %g\n', mtc.sefoffset);
   fprintf(fo,'#datatype: %d\n', mtc.datatype);
   fprintf(fo,'#Each row i contains %d values for node i\n', mtc.numvol);
   Opt1D.Space = 't';
   Opt1D.Fast = 'y';
   Opt1D.verbose = 0;
   Opt1D.OverWrite = 'a';
   [err, UsedName] = wryte3(mtc.data, o1d, Opt1D);
   [h,m,s] = sectotime(etime(clock, to));
   fprintf(2,'Writing took %g:%g:%g.\n', h,m,s);
end
