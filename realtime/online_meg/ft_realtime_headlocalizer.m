function ft_realtime_headlocalizer(cfg)

% FT_REALTIME_HEADLOCALIZER is a realtime application for online
% visualization of the head localization coils in a CTF275 system.
%
% Use as
%   ft_realtime_headlocalizer(cfg)
% with the following configuration options
%   cfg.template        = string, name of the original data set to be used as template (default = [])
%   cfg.blocksize       = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.bufferdata      = whether to start on the 'first or 'last' data that is available (default = 'last')
%   cfg.accuracy_green  = distance from fiducial coordinate; green when within limits (default =  0.15 cm)
%   cfg.accuracy_orange = orange when within limits, red when out (default =  0.3 cm)
%
% The source of the data is configured as
%   cfg.dataset       = string, default is 'buffer://odin:1972'
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic


% Copyright (C) 2008-2012, Robert Oostenveld, Arjen Stolk, Ana Todorovic, Stefan Klanke
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% defaults
ft_defaults
cfg.template        = ft_getopt(cfg, 'template',         []); % template dataset containing the references
cfg.blocksize       = ft_getopt(cfg, 'blocksize',         1); % in seconds
cfg.bufferdata      = ft_getopt(cfg, 'bufferdata',   'last'); % first (replay) or last (real-time)
cfg.dataset         = ft_getopt(cfg, 'dataset', 'buffer://odin:1972'); % location of the buffer/dataset
cfg.accuracy_green  = ft_getopt(cfg, 'accuracy_green',  .15); % green when within this distance from reference
cfg.accuracy_orange = ft_getopt(cfg, 'accuracy_orange',  .3); % orange when within this distance from reference

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

% read the template coil positions
if ~isempty(cfg.template)
    template = ft_read_headshape(cfg.template, 'coordinates', 'dewar');
    reference(1,:) = [template.fid.pnt(1,1), template.fid.pnt(1,2), template.fid.pnt(1,3)]; % nasion
    reference(2,:) = [template.fid.pnt(2,1), template.fid.pnt(2,2), template.fid.pnt(2,3)];
    reference(3,:) = [template.fid.pnt(3,1), template.fid.pnt(3,2), template.fid.pnt(3,3)];
else
    reference = [];
end

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true);

% MEG sensors in CTF dewar space
sens = headcoordinates2ctfdewar(hdr.orig.hc.dewar(:,1)', hdr.orig.hc.dewar(:,2)', hdr.orig.hc.dewar(:,3)', hdr.grad);

% define a subset of channels for reading, only "headloc" type channels are relevant
if strcmp(cfg.dataset, 'buffer://odin:1972');
    chanindx = 1:9; % odin buffer specific
else
    [~, chanindx] = match_str('headloc', ft_chantype(hdr));
end

if isempty(chanindx)
    error('the data does not seem to have headlocalization channels');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
prevSample  = 0;
count       = 0;

% initiate figure
hMainFig = figure;
set(hMainFig, 'UserData', 1); % initialize the flag variable
set(hMainFig, 'KeyPressFcn', {@key_sub});

% insert gui variables
info           = [];
info.cfg       = cfg;
info.reference = reference;
info.sens      = sens;
guidata(hMainFig, info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (get(hMainFig, 'UserData') == 1) % while the flag is one, the loop continues
    
    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.headerfile, 'cache', true);
    
    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);
    
    if newsamples>=blocksize
        
        if strcmp(cfg.bufferdata, 'last')
            begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
            endsample  = hdr.nSamples*hdr.nTrials;
        elseif strcmp(cfg.bufferdata, 'first')
            begsample  = prevSample + 1;
            endsample  = prevSample + blocksize ;
        else
            error('unsupported value for cfg.bufferdata');
        end
        
        % remember up to where the data was read
        prevSample  = endsample;
        count       = count + 1;
        fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
        
        % read data segment from buffer
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here onward it is specific to the headlocalization in the CTF system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % put the data in a fieldtrip-like raw structure
        data.trial{1} = double(dat);
        data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
        data.label    = hdr.label(chanindx);
        data.hdr      = hdr;
        data.fsample  = hdr.Fs;
        clear dat;
        
        % assign the channels to the resp. coil coordinates
        [~, x1] = match_str('HLC0011', data.label);
        [~, y1] = match_str('HLC0012', data.label);
        [~, z1] = match_str('HLC0013', data.label);
        [~, x2] = match_str('HLC0021', data.label);
        [~, y2] = match_str('HLC0022', data.label);
        [~, z2] = match_str('HLC0023', data.label);
        [~, x3] = match_str('HLC0031', data.label);
        [~, y3] = match_str('HLC0032', data.label);
        [~, z3] = match_str('HLC0033', data.label);
        
        % get gui variables
        info = guidata(hMainFig);
        
        % convert from meter to cm and assign to the resp. coil
        info.coil1 = data.trial{1}([x1 y1 z1],:) * 100;
        info.coil2 = data.trial{1}([x2 y2 z2],:) * 100;
        info.coil3 = data.trial{1}([x3 y3 z3],:) * 100;      
        
        % store gui variables
        guidata(hMainFig, info);
        
        % DRAW LEFT PANEL - BACK VIEW
        a = subplot(1,2,1);
        h = get(a, 'children');
        hold on;
        
        if ~isempty(h)
            % done on every iteration
            delete(h);
        end
        
        % draw the color-coded head and distances from the references
        draw_sub(hMainFig);
        
        % show the update key
        str = sprintf('Press "u" to update the reference coordinates \n Press "s" to toggle on/off sensors/dewar display \n Press "q" to stop the application \n');
        title(str);
        
        if isempty(h)
            % only done once
            grid on
            xlabel('x (cm)');
            ylabel('y (cm)');
            zlabel('z (cm)');
            set(gca, 'xtick', -20:2:20)
            set(gca, 'ytick', -20:2:20)
            set(gca, 'ztick', -50:2:0) % note the different scaling
            view(-45, 90)
            % axis square
            axis vis3d
            axis manual
        end
        
        % DRAW RIGHT PANEL - TOP VIEW
        b = subplot(1,2,2);
        i = get(b, 'children');
        hold on;
        
        if ~isempty(i)
            % done on every iteration
            delete(i);
        end
        
        % draw the color-coded head and distances from the references
        draw_sub(hMainFig);
        
        % show current timesample
        str = sprintf('Runtime = %d s\n', round(mean(data.time{1})));
        clear data;
        title(str);
        fprintf(str);
        
        if isempty(i)
            % only done once
            grid on
            xlabel('x (cm)');
            ylabel('y (cm)');
            zlabel('z (cm)');
            set(gca, 'xtick', -20:2:20)
            set(gca, 'ytick', -20:2:20)
            set(gca, 'ztick', -50:2:0) % note the different scaling
            view(-45, 0)
            % axis square
            axis vis3d
            axis square
        end
        
        % force Matlab to update the figure
        drawnow
        
    end % if enough new samples
end % while true
close(hMainFig); % close the figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which transform the positions from headcoordinate to dewar space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_grad = headcoordinates2ctfdewar(nas, lpa, rpa, grad) % based on headcoordinates.m

% compute the origin and direction of the coordinate axes in MRI coordinates
origin = (lpa+rpa)/2;
dirx = nas-origin;
dirx = dirx/norm(dirx);
dirz = cross(dirx,lpa-rpa);
dirz = dirz/norm(dirz);
diry = cross(dirz,dirx);

% compute the rotation matrix
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);

% compute the translation matrix
tra = eye(4);
tra(1:4,4)   = [-origin(:); 1];

% compute the full homogeneous transformation matrix from these two
h = rot * tra;

% we do not want to plot the REF sensors
chansel = match_str(grad.chantype,'meggrad');
grad.chanpos  = grad.chanpos(chansel,:);
grad.chanori  = grad.chanori(chansel,:);
grad.chantype = grad.chantype(chansel,:);
grad.label    = grad.label(chansel,:);

% apply the inverse transformation matrix on data
t_grad = grad;
t_grad.chanpos = warp_apply(inv(h), grad.chanpos, 'homogenous');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does the timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time] = offset2time(offset, fsample, nsamples)

offset   = double(offset);
nsamples = double(nsamples);
time = (offset + (0:(nsamples-1)))/fsample;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which computes the circumcenter(x,y,z) of the 3D triangle (3 coils)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cc] = circumcenter(coil1, coil2, coil3)

% use coordinates relative to point `a' of the triangle
xba = coil2(1,end) - coil1(1,end);
yba = coil2(2,end) - coil1(2,end);
zba = coil2(3,end) - coil1(3,end);
xca = coil3(1,end) - coil1(1,end);
yca = coil3(2,end) - coil1(2,end);
zca = coil3(3,end) - coil1(3,end);

% squares of lengths of the edges incident to `a'
balength = xba * xba + yba * yba + zba * zba;
calength = xca * xca + yca * yca + zca * zca;

% cross product of these edges
xcrossbc = yba * zca - yca * zba;
ycrossbc = zba * xca - zca * xba;
zcrossbc = xba * yca - xca * yba;

% calculate the denominator of the formulae
denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc);

% calculate offset (from `a') of circumcenter
xcirca = ((balength * yca - calength * yba) * zcrossbc - (balength * zca - calength * zba) * ycrossbc) * denominator;
ycirca = ((balength * zca - calength * zba) * xcrossbc - (balength * xca - calength * xba) * zcrossbc) * denominator;
zcirca = ((balength * xca - calength * xba) * ycrossbc - (balength * yca - calength * yba) * xcrossbc) * denominator;

cc(1) = xcirca + coil1(1,end);
cc(2) = ycirca + coil1(2,end);
cc(3) = zcirca + coil1(3,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which draws the color-coded head and distances to the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_sub(handle)

% get the info
info      = guidata(handle);

if strcmp(get(handle,'Tag'), 'On')
    % plot the sensors
    hold on; ft_plot_sens(info.sens);
end

% plot the three reference fiducial positions
if ~isempty(info.reference)
    plot3(info.reference(1,1), info.reference(1,2), info.reference(1,3), 'k^', 'MarkerSize', 27, 'LineWidth', 2);
    plot3(info.reference(2,1), info.reference(2,2), info.reference(2,3), 'ko', 'MarkerSize', 27, 'LineWidth', 2);
    plot3(info.reference(3,1), info.reference(3,2), info.reference(3,3), 'ko', 'MarkerSize', 27, 'LineWidth', 2);
    
    text(-8,8, info.reference(2,3), 'Left', 'FontSize', 15);
    text(6,-6, info.reference(3,3), 'Right', 'FontSize', 15);
end

% plot the coil position traces
plot3(info.coil1(1,:), info.coil1(2,:), info.coil1(3,:), 'k^', 'LineWidth', 1,'MarkerSize', 3) % nasion
plot3(info.coil2(1,:), info.coil2(2,:), info.coil2(3,:), 'ko', 'LineWidth', 1,'MarkerSize', 3)
plot3(info.coil3(1,:), info.coil3(2,:), info.coil3(3,:), 'ko', 'LineWidth', 1,'MarkerSize', 3)

% draw nasion position
if ~isempty(info.reference)
    if abs(info.reference(1,1))-info.cfg.accuracy_green < abs(info.coil1(1,end)) && abs(info.coil1(1,end)) < abs(info.reference(1,1))+info.cfg.accuracy_green ...
            && abs(info.reference(1,2))-info.cfg.accuracy_green < abs(info.coil1(2,end)) && abs(info.coil1(2,end)) < abs(info.reference(1,2))+info.cfg.accuracy_green ...
            && abs(info.reference(1,3))-info.cfg.accuracy_green < abs(info.coil1(3,end)) && abs(info.coil1(3,end)) < abs(info.reference(1,3))+info.cfg.accuracy_green
        plot3(info.coil1(1,end),info.coil1(2,end),info.coil1(3,end),'g^', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        head1 = true;
    elseif abs(info.reference(1,1))-info.cfg.accuracy_orange < abs(info.coil1(1,end)) && abs(info.coil1(1,end)) < abs(info.reference(1,1))+info.cfg.accuracy_orange ...
            && abs(info.reference(1,2))-info.cfg.accuracy_orange < abs(info.coil1(2,end)) && abs(info.coil1(2,end)) < abs(info.reference(1,2))+info.cfg.accuracy_orange ...
            && abs(info.reference(1,3))-info.cfg.accuracy_orange < abs(info.coil1(3,end)) && abs(info.coil1(3,end)) < abs(info.reference(1,3))+info.cfg.accuracy_orange
        plot3(info.coil1(1,end),info.coil1(2,end),info.coil1(3,end),'y^', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        head1 = false;
    else % when not in correct position
        plot3(info.coil1(1,end),info.coil1(2,end), info.coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        head1 = false;
    end
else
    plot3(info.coil1(1,end),info.coil1(2,end), info.coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
    head1 = false;
end

% draw left ear position
if ~isempty(info.reference)
    if abs(info.reference(2,1))-info.cfg.accuracy_green < abs(info.coil2(1,end)) && abs(info.coil2(1,end)) < abs(info.reference(2,1))+info.cfg.accuracy_green ...
            && abs(info.reference(2,2))-info.cfg.accuracy_green < abs(info.coil2(2,end)) && abs(info.coil2(2,end)) < abs(info.reference(2,2))+info.cfg.accuracy_green ...
            && abs(info.reference(2,3))-info.cfg.accuracy_green < abs(info.coil2(3,end)) && abs(info.coil2(3,end)) < abs(info.reference(2,3))+info.cfg.accuracy_green
        plot3(info.coil2(1,end),info.coil2(2,end),info.coil2(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        head2 = true;
    elseif abs(info.reference(2,1))-info.cfg.accuracy_orange < abs(info.coil2(1,end)) && abs(info.coil2(1,end)) < abs(info.reference(2,1))+info.cfg.accuracy_orange ...
            && abs(info.reference(2,2))-info.cfg.accuracy_orange < abs(info.coil2(2,end)) && abs(info.coil2(2,end)) < abs(info.reference(2,2))+info.cfg.accuracy_orange ...
            && abs(info.reference(2,3))-info.cfg.accuracy_orange < abs(info.coil2(3,end)) && abs(info.coil2(3,end)) < abs(info.reference(2,3))+info.cfg.accuracy_orange
        plot3(info.coil2(1,end),info.coil2(2,end),info.coil2(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        head2 = false;
    else % when not in correct position
        plot3(info.coil2(1,end),info.coil2(2,end), info.coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        head2 = false;
    end
else
    plot3(info.coil2(1,end),info.coil2(2,end), info.coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
    head2 = false;
end

% draw right ear position
if ~isempty(info.reference)
    if abs(info.reference(3,1))-info.cfg.accuracy_green < abs(info.coil3(1,end)) && abs(info.coil3(1,end)) < abs(info.reference(3,1))+info.cfg.accuracy_green  ...
            && abs(info.reference(3,2))-info.cfg.accuracy_green  < abs(info.coil3(2,end)) && abs(info.coil3(2,end)) < abs(info.reference(3,2))+info.cfg.accuracy_green  ...
            && abs(info.reference(3,3))-info.cfg.accuracy_green  < abs(info.coil3(3,end)) && abs(info.coil3(3,end)) < abs(info.reference(3,3))+info.cfg.accuracy_green
        plot3(info.coil3(1,end),info.coil3(2,end),info.coil3(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        head3 = true;
    elseif abs(info.reference(3,1))-info.cfg.accuracy_orange < abs(info.coil3(1,end)) && abs(info.coil3(1,end)) < abs(info.reference(3,1))+info.cfg.accuracy_orange ...
            && abs(info.reference(3,2))-info.cfg.accuracy_orange < abs(info.coil3(2,end)) && abs(info.coil3(2,end)) < abs(info.reference(3,2))+info.cfg.accuracy_orange ...
            && abs(info.reference(3,3))-info.cfg.accuracy_orange < abs(info.coil3(3,end)) && abs(info.coil3(3,end)) < abs(info.reference(3,3))+info.cfg.accuracy_orange
        plot3(info.coil3(1,end),info.coil3(2,end),info.coil3(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        head3 = false;
    else % when not in correct position
        plot3(info.coil3(1,end),info.coil3(2,end), info.coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        head3 = false;
    end
else
    plot3(info.coil3(1,end),info.coil3(2,end), info.coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
    head3 = false;
end

% draw 3d head
cc = circumcenter(info.coil1,info.coil2,info.coil3);
x_radius = sqrt((info.coil2(1,end) - cc(1))^2 + (info.coil2(2,end) - cc(2))^2);
y_radius = sqrt((info.coil3(1,end) - cc(1))^2 + (info.coil3(2,end) - cc(2))^2);
[xe, ye, ze] = ellipsoid(cc(1),cc(2),cc(3),x_radius,y_radius,11);
hh = surfl(xe, ye, ze);
shading interp
if head1 == true && head2 == true && head3 == true
    colormap cool
else
    colormap hot
end
alpha(.15)

% put the info back
guidata(handle, info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata)

% get the info
info = guidata(handle);

switch eventdata.Key
    case 'u'
        % update the info.reference positions
        fprintf('updating reference coordinates \n')
        info.reference(1,:) = info.coil1(:,end); % nasion
        info.reference(2,:) = info.coil2(:,end); % right ear
        info.reference(3,:) = info.coil3(:,end); % left ear
    case 's'
        % display the sensors/dewar
        if isempty(get(1,'Tag'))
            fprintf('displaying sensors/dewar \n')
            set(handle,'Tag', 'On'); % toggle on
        elseif strcmp(get(1,'Tag'), 'On')
            set(handle,'Tag', ''); % toggle off
        end
    case 'q'
        % stop the application
        fprintf('stopping the application \n')
        set(handle, 'UserData', 0); % stop the while loop
    otherwise
        fprintf('no command executed \n')
end

% put the info back
guidata(handle, info);