function replay_dicoms(src, dest)
% function replay_dicoms(src, dest)
% src = source directory (where the dicom files are)
% dest = destination directory (mrprot.txt + pixeldata)

FNs = dir([src '/*.ima']);
N = numel(FNs);

if N<1
   error('No DICOM files found');
end

fprintf(1,'Found %i DICOM files\n', N);

I1 = dicominfo([src '/' FNs(1).name]);
if ~isfield(I1, 'Private_0029_1020')
	error('Private field 0029:1020 not found - make sure this comes from a Siemens scanner');
end
prot = char(I1.Private_0029_1020)';

begpos = strfind(prot, '### ASCCONV BEGIN ###');
endpos = strfind(prot, '### ASCCONV END ###');

if numel(begpos)~=1 || numel(endpos)~=1
   error('ASCII protocol information not found');
end

mrprot = prot(begpos+21:endpos-1);

F = fopen([dest '/mrprot.txt'],'w');
fwrite(F, mrprot, 'char');
fclose(F);
T0 = tic;

for n=1:N
	D = dicomread([src '/' FNs(n).name])';
	imagesc(D);
	drawnow;
	
	pause(2.0*n - toc(T0));
	dname = sprintf('%s/%06i.PixelData', dest, n);
	F = fopen(dname, 'wb');
	fwrite(F, D, 'int16');
	fclose(F);
	fprintf(1,'Wrote out scan %i\n', n);
end
   

