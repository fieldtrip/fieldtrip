% CPD 2D image registration example.
% The idea is to extract feature points from both images beforehand.
% Register the extracted point sets. And then apply the estimated
% transformation to the images.

clear all; close all; clc;
load cpd_data_geometry.mat; % load two images

% Extract feature points. Here we use a simple edge detector.
% Most likely you'll have to use a more advanced feature detector here.
[j,i]=find(edge(im,'canny'));
X=[i j]; % first point set

[j,i]=find(edge(refim,'canny'));
Y=[i j]; % second point set


% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=1;          % show every iteration

% registering Y to X
Transform=cpd_register(X,Y,opt);

% Now lets apply the found transformation to the image
% Create a dense grid of points
[M,N]=size(refim);
[x,y]=meshgrid(1:N,1:M); 
grid=[x(:) y(:)];

% Transform the grid according to the estimated transformation
T=cpd_transform(grid, Transform);

% Interpolate the image
Tx=reshape(T(:,1),[M N]);
Ty=reshape(T(:,2),[M N]);
result=interp2(im,Tx,Ty);   % evaluate image Y at the new gird (through interpolation)

% Show the result
figure,imshow(im); title('Source image');
figure,imshow(refim); title('Reference image');
figure,imshow(result); title('Result image');
