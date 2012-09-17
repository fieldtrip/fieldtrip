function test_triangulate_seg

% TEST: test_triangulate_seg
% TEST: triangulate_seg

[ftpath,n,e] = fileparts(which('ft_defaults'));
pwdir  = pwd;
cd(ftpath);
cd('private');

% create a segmented volume containing 2 blobs with a hole in it
seg = false(100,100,100);
seg(21:80,21:80,21:80) = true;
seg(3:15,3:15,3:15)    = true;
seg(31:50,31:50,31:50) = false;
[pnt, tri] = triangulate_seg(seg, 500); 

assert(all(round(mean(pnt))==50));

cd(pwdir);
