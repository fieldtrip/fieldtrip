function dat = readmptiff(fname)
%
% vol=readmptiff(fname)
%
% load a volume from a multi-page TIFF file
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%      fname: input file name
%
% output:
%      dat: output, data read from the TIFF file
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

info = imfinfo(fname);
nslice = length(info);
if (nslice <= 0)
    error('no data found in the tiff');
end

for i = 1:nslice
    dat(:, :, i) = imread(fname, i);
end
