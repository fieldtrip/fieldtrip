function test_bug1976

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1976 ft_sourceanalysis test_ft_sourceanalysis

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% clear all
% close all
% addpath('/users/MEG-Data/misc/pawel/tmp/fieldtrip-20130205/')
% ft_defaults

%% simple head model
N_phi = 100;
N_theta = 100;
R = 10;
[phi theta] = meshgrid(linspace(0,2*pi,N_phi), linspace(-pi/2,pi/2,N_theta));
index = theta == pi/2;
index(find(index,1)) = false;
theta(index) = [];
phi(index) = [];
index = theta == -pi/2;
index(find(index,1)) = false;
theta(index) = [];
phi(index) = [];
pnt = [R*cos(theta(:)).*cos(phi(:)) R*cos(theta(:)).*sin(phi(:)) R*sin(theta(:))];
vol = ft_headmodel_singlesphere(pnt);

%% a few sensors in old FieldTrip style for magnetometer
N_phi = 15;
N_theta = 15;
R = 12;
[phi theta] = meshgrid(linspace(0,2*pi,N_phi), linspace(0,pi/2,N_theta));
index = theta == pi/2;
index(find(index,1)) = false;
theta(index) = [];
phi(index) = [];
sens.pnt = [R*cos(theta(:)).*cos(phi(:)) R*cos(theta(:)).*sin(phi(:)) R*sin(theta(:))];
sens.ori = sens.pnt./repmat(sqrt(sum(sens.pnt.^2,2)),1,3);
sens.label    = cell(size(sens.pnt,1),1);
for ii = 1:length(sens.label)
    sens.label{ii} = num2str(ii);
end
[vol, sens] = ft_prepare_vol_sens(vol, sens);
%% siple grid
N_phi = 5;
N_theta = 5;
N_R = 5;
[phi theta R] = meshgrid(linspace(0,2*pi,N_phi), linspace(-pi/2,pi/2,N_theta), linspace(1,9,N_R));
for r = linspace(1,9,N_R)
    index = theta == pi/2 & R == r;
    index(find(index,1)) = false;
    theta(index) = [];
    phi(index) = [];
    R(index) = [];
    index = theta == -pi/2 & R == r;
    index(find(index,1)) = false;
    theta(index) = [];
    phi(index) = [];
    R(index) = [];
end
pos = [R(:).*cos(theta(:)).*cos(phi(:)) R(:).*cos(theta(:)).*sin(phi(:)) R(:).*sin(theta(:))];

%% leadfield
lf = ft_compute_leadfield(pos, sens, vol, 'normalizeparam', 0);

%% forward
j = rand(3,size(pos,1));
v_r = (pos./repmat(sqrt(sum(pos.^2,2)),1,3)).';
j = j - v_r.*repmat(sum(j.*v_r),3,1);
signal = lf * j(:);

data = [];
data.time = 0;
data.label = sens.label;
data.avg = signal;
data.dimord = 'chan_time';
%% inverse
cfg = [];
cfg.method = 'mne';
cfg.grid.pos = pos;
cfg.grid.inside = 1:size(cfg.grid.pos,1);
cfg.grid.outside = [];
cfg.grid.leadfield = cell(1,size(pos,1));
for ii = 1:length(cfg.grid.leadfield)
    cfg.grid.leadfield{ii} = lf(:,(1:3) + (ii - 1) * 3);
end
cfg.mne.snr = 1;
cfg.vol = vol;
cfg.grad = sens;
source = ft_sourceanalysis(cfg, data);

%% test plot
plot3(pnt(:,1),pnt(:,2),pnt(:,3),'y.'); hold on
quiver3(sens.chanpos(:,1),sens.chanpos(:,2),sens.chanpos(:,3),sens.chanori(:,1),sens.chanori(:,2),sens.chanori(:,3))
quiver3(pos(:,1),pos(:,2),pos(:,3),j(1,:)',j(2,:)',j(3,:)')

%%
% Note: I (Johanna) tested this also by creating a temporary timelock ctf151 with
% only 1 trial and 1 time point and running test_ft_sourceanalysis with
% this dataset (see also commented out section of ref_datasets) but as this
% does not make sense for the freq cases, this would not normally run as
% part of the test_ft_sourceanalysis function.

