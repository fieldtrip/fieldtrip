function filename = duneuro_installbrainstormapp(where)

if nargin<1
  where = tempdir;
end 

url     = 'https://neuroimage.usc.edu/bst/getupdate.php?d=bst_duneuro.zip';
zipFile = fullfile(where, 'bst_duneuro.zip');
outDir  = fullfile(where);  % change this to your desired target folder

% Download
websave(zipFile, url);

% Unzip
unzip(zipFile, outDir);

% (Optional) clean up
delete(zipFile);

where = fullfile(where, 'bst_duneuro', 'bin');

str = computer('arch');
switch str
  case 'win64'
    filename = 'bst_duneuro_meeg_win64.exe';
  case 'glnxa64'
    filename = 'bst_duneuro_meeg_linux64.app';
  case 'maci64'
    filename = 'bst_duneuro_meeg_mac64.app';
  case 'maca64'
    error('no executable for this platform');
  otherwise
    error('no executable for this platform');
end

filename = fullfile(where, filename);
