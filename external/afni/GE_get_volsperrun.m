function [ volsperrun ] = GE_get_volsperrun(inDir)
%
% [ volsperrun ] = GE_get_volsperrun(inDir)
%
% A utility to find the number of volumes per run in an ADW scan.
% inDir is the starting directory for the I.* files.
%
% Souheil J. Inati  
% Dartmouth College
% July 2001
% souheil.inati@dartmouth.edu
%

% Get series ID number
ser_idstr = inDir(end);
ser_idnum = str2num(ser_idstr);

% Strip tailing directory
GE_data_root = inDir(1:max(findstr('/',inDir))-1);

% List directories and extract those corresponding to the desired series
d = dir(GE_data_root);

dir_idx = find(cat(d.isdir));
ndir = length(dir_idx);

dirmat = char(d.name);

good_idx = find((dirmat(:,3) == ser_idstr) & isspace(dirmat(:,4)));

ser_dirnames = char(d(good_idx).name);

ndir = size(ser_dirnames,1);

byte_align = 0;

% Open the first file and get the number of slices
firstfile = fullfile(GE_data_root, ser_dirnames(1,:), 'I.001');
[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(firstfile);
nslices = im_hdr.slquant;

nimg = 0;

for idir = 1:ndir
  curdir = fullfile(GE_data_root, ser_dirnames(idir,:));
  subd = dir(curdir);
  
  dirmat = char(subd.name);
  good_idx = strmatch('I.', dirmat);
  nfiles = length(good_idx);
  
  disp(sprintf('Found %d image files in directory %s', nfiles, curdir))
  
  flist = sort(cellstr(char(subd(good_idx).name)));

  % Loop through and read the run number from the header of each file
  for ifile = 1:nfiles
    imageFile = fullfile(curdir, flist{ifile});
    nimg = nimg+1;
    fid = fopen(imageFile,'r','b');
    fseek(fid,4936,'bof');
    runnum(nimg) =  fread(fid,1,'float32'); % im_hdr.user17
    imagename{nimg} = imageFile;
    fclose(fid);
  end
end

% Output the results to the screen
[B,I,J] = unique(runnum);

imperrun = [I(1) diff(I)];

volsperrun = imperrun/nslices;

%numrun = length(B);
%fprintf('There are %d functionals runs.\n',numrun)

%for i = 1:numrun
%  fprintf('run %d consists of %d images = %d volumes\n',B(i),imperrun(i),imperrun(i)/nslices)
%end

return



    %    im_rawrunnum(nimg) = im_hdr.rawrunnum;
    %    im_datetime(nimg) = im_hdr.im_datetime;
    %    im_actual_dt(nimg) = im_hdr.im_actual_dt;
    %    im_lastmod(nimg) = im_hdr.im_lastmod;
    %    im_checksum(nimg) = im_hdr.im_checksum;
    %    im_scanactno(nimg) = im_hdr.scanactno;


