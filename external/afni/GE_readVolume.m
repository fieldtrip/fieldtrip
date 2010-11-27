function [imageVol, lastfile] = GE_readVolume(startDir, passnum, volSize, depth, im_offset)
%
%GE_readVolume
% 
% [imageVol, lastfile] = GE_readVolume(startDir, passnum, [nX nY nZ], depth, im_offset)
%
% reads the volume for passnum from the series which is stored
% starting in startDir and returns the name of the last file read
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

% initialize some variables
nX = volSize(1);
nY = volSize(2);
nZ = volSize(3);
sliceSize = nX*nY;
imageVol = zeros(nX, nY, nZ);
[path_stem, path_start] = fileparts(startDir);

for i = 1:nZ
        % Make the filename
        filenum = (passnum-1)*nZ + i;   % file no.
        filep = ceil(filenum/999) - 1;
        filen = filenum - 999*filep;
        path_num = str2num(path_start) + 20*filep;
        path_now = sprintf('00%d',path_num);
        path_now = path_now(length(path_now)-2:length(path_now));
        path = fullfile(path_stem, path_now);
	stub = sprintf('00%d',filen);
	stub = stub(length(stub)-2:length(stub));
        imageFile = fullfile(path,['I.' stub]);

        % Open the file
	[fid,message] = fopen(imageFile,'r','b');
	if (fid == -1)
          fprintf('Cannot Open %s.\n',imageFile);
          passnum
	  break
	end

	% Skip to the data
	fseek(fid,im_offset,-1);
	% Read the slice
	buffer = fread(fid,sliceSize,sprintf('int%d',depth));
	% append the slice to the imageSet
	imageVol(:,:,i) = reshape(buffer, nX, nY);

        % Close the file
	status = fclose(fid);

        % Set the lastfile
        lastfile = imageFile;
end

return
