function test_triangulate_seg

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY triangulate_seg
% DATA no

% since the function to test is in a private directory, we explicitly have to cd into that directory
[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'private'));

% create a segmented volume containing 2 blobs with a hole in it
seg = false(100,100,100);
seg(21:80,21:80,21:80) = true;
seg(3:15,3:15,3:15)    = true;
seg(31:50,31:50,31:50) = false;
[pnt, tri] = triangulate_seg(seg, 500); 

assert(all(round(mean(pnt))==50));

