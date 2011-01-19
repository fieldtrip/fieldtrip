clear;

%% Load the reciprocity data
infile = '../data/SPM/recipdata-gray.h5';
[data,compress,gridlocs,node_sizes,voxel_sizes] = ...
    hh_read_recipdata(infile);
decompress = gridlocs(:,1) + gridlocs(:,2) * node_sizes(1) + ...    
    gridlocs(:,3) * node_sizes(2) * node_sizes(1) + 1;
dipoles = hdf5read('test/dipoles.h5','/model/dipoles');

%% Obtain the lead field G for all dipoles
NDIPOLE = size(dipoles,1);
N = size(data,1)

%% Generate the test data 
% $$$ Generate the dipole sources
NSAMPLE = 200;
T = 1e-3;
A1 = 4e-6;
A2 = 4e-6;
A3 = 1e-6;
A4 = 1e-6;
t = 1:NSAMPLE;
f1 = 15;
f2 = 18;
f3 = 15;
f4 = 15;
phi1 = 0;
phi2 = pi / 3;
phi3 = 0;
phi4 = 2 * pi / 3;

x1 = A1 * sin(2 * f1 * pi * t * T + phi1);
x2 = A2 * sin(2 * f2 * pi * t * T + phi2);
x3 = A3 * sin(2 * f3 * pi * t * T + phi3);
x4 = A4 * sin(2 * f4 * pi * t * T + phi4);

% $$$ Compute the lead field matrices for two dipoles
pos1 = [119,27,156];
pos2 = [119,176,179];
pos3 = [119,105,222];
pos4 = [72,105,222];

eta1 = [1,0,0]';
eta2 = [-1,0,0]';
eta3 = [0,0,1]';
eta4 = [0,0,1]';

L1 = hh_leadfield(data,compress,node_sizes,pos1(1),pos1(2), ...
                  pos1(3));
L2 = hh_leadfield(data,compress,node_sizes,pos2(1),pos2(2), ...
                  pos2(3));
L3 = hh_leadfield(data,compress,node_sizes,pos3(1),pos3(2), ...
                  pos3(3));
L4 = hh_leadfield(data,compress,node_sizes,pos4(1),pos4(2), ...
                  pos4(3));

g1 = L1 * eta1;
g2 = L2 * eta2;
g3 = L3 * eta3;
g4 = L4 * eta4;

% $$$ Obtain the simulated brain wave
% $$$ m = g1 * x1 + g2 * x2 + g3 * x3 + g4 * x4;
m = g1 * x1 + g2 * x2 + g3 * x3;
% $$$ m = g3 * x3 + g4 * x4;
snr = 16;
y = awgn(m,snr,'measured','db');
plot(y');

% Compute the covariance matrix
R = y * y';
N = length(R);
[U,S,V] = svd(R);
gamma = 0.001;
lambda_max = max(eig(R));
Sigma_inv = inv(R + eye(N,N) * gamma * lambda_max);

%% Obtain the spectral covariance matrix

%% Compute the Backus-Gilbert covariance matrix
% $$$ G = hh_region_weight(data,compress,node_sizes,dipoles,[]);
newrgn = hh_sphere_region(dipoles,10,pos3);
G = hh_region_weight(data,compress,node_sizes,newrgn,[]);
Rw = G * G';
scale = norm(R,'fro') / norm(Rw,'fro')
Rw = scale * Rw;
lambda = 1;
Rnew = R + lambda * Rw;
gamma = 0.001;
lambda_max = max(eig(Rnew));
Rnew_inv = inv(Rnew + eye(N,N) * gamma * lambda_max);

%% Test distributed source approach
newrgn = hh_sphere_region(dipoles,30,pos4);
G = hh_region_weight(data,compress,node_sizes,newrgn,[]);

%% Obtain the fMRI-guided covariance matrix

%% Obtain the tomographic map of LCMV
NSOURCE = 4;
Nk = [];
% $$$ cost = hh_music(data,compress,decompress,node_sizes,R,NSOURCE,Nk);
cost_lcmv = hh_lcmv_map(data,compress,decompress,node_sizes, R, ...
                        Sigma_inv);

%% V-BKB
cost_vbkb = hh_vbkb_map(data,compress,decompress,node_sizes, R, ...
                        Sigma_inv);

%% Compute the tomographic map of SAM
cost_sam = hh_sam_map(data,compress,decompress,node_sizes, R, ...
                      Sigma_inv);

%% fMRI guided beamformer
cost_fmri = hh_sam_map(data,compress,decompress,node_sizes, Rnew, ...
                       Rnew_inv);

%% Compute the tomographic map of SAM using SVD approach
cost_samsvd = hh_sam_mapsvd(data,compress,decompress,node_sizes, R, ...
                            Sigma_inv);

%% Obtaine the tomographic map of DICS

%% Obtain the tomographic map of V-WNMVB
cost_vwnmvb = hh_vwnmvb_map(data,compress,decompress,node_sizes, R, ...
                            Sigma_inv);

%% Obtain the tomographic map of V-WNMVB using SVD approach
cost_vwnmvbsvd = hh_vwnmvb_mapsvd(data,compress,decompress,node_sizes, R, ...
                                  Sigma_inv);

%% Compute the tomographic map of MUSIC
cost_music = hh_music(data,compress,decompress,node_sizes,R,4,[]);

%% Save the solution to the output file to create the 3D volume
%% image
btype = 'fmri';
switch btype
  case 'sam'
    cost = cost_sam;
    tail = 'sam';
  case 'samsvd'
    cost = cost_sam;
    tail = 'samsvd';
  case 'lcmv'
    cost = cost_lcmv;
    tail = 'lcmv';
  case 'vwnmvb'
    cost = cost_vwnmvb;
    tail = 'vwnmvb';
  case 'vwnmvbsvd'
    cost = cost_vwnmvbsvd;
    tail = 'vwnmvbsvd';
  case 'vbkb'
    cost = cost_vbkb;
    tail = 'vbkb';
  case 'music'
    cost = cost_music;
    tail = 'music';    
  case 'fmri'
    cost = cost_fmri;
    tail = 'fMRI';        
  otherwise
end

outfile = 'tmp/fdm_paper.h5';
newcompress = zeros(length(compress),1);
newcompress(decompress) = 1:length(decompress);
hdf5write(outfile,'/data/info',['This file containts the output ' ...
                    'solution for the FDN paper']);
hdf5write(outfile,'/data/sol',cost,'WriteMode', 'append');
hdf5write(outfile,'/data/compress',int32(newcompress - 1),'WriteMode', 'append');
hdf5write(outfile,'/data/node_sizes',int32(node_sizes),'WriteMode', ...
          'append');

% $$$ Convert the solution to the ANALYZE volume image
str = sprintf(['/home/hungptit/working/fdm/bin/solimg --datafile ' ...
               'tmp/fdm_paper.h5 --alynzfile tmp/fdn_%s'],tail);
% $$$ unix(['/home/hungptit/working/fdm/bin/solimg --datafile tmp/fdm_paper.h5 ' ...
% $$$       '--alynzfile tmp/fdn']);
unix(str);
hh_create_colormap('tmp/fdn',0.1,0);

%% View solution with MATLAB
% Decompress the solution
X = zeros(node_sizes(1),node_sizes(2),node_sizes(3));
X(decompress) = cost;
% Display the output map
imagesc(squeeze(X(:,:,222)));