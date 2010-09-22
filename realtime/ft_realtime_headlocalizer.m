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
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% To stop the realtime function, you have to press Ctrl-C

% Copyright (C) 2008-2010, Robert Oostenveld, Arjen Stolk, Ana Todorovic, Stefan Klanke
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information


if ~isfield(cfg, 'template'),       cfg.template = [];        end
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last
% distance from fiducial coordinate that makes the localizer turn green or orange
if ~isfield(cfg, 'accuracy_green'),       cfg.accuracy_green = .15;        end
if ~isfield(cfg, 'accuracy_orange'),      cfg.accuracy_orange = .3;        end

% translate dataset into datafile+headerfile
cfg = checkconfig(cfg, 'dataset2files', 'yes');

% read the template coil positions
if ~isempty(cfg.template)
    template = read_headshape(cfg.template, 'coordinates', 'dewar');
else
    template = [];
end

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading, only "headloc" type channels are relevant
if strcmp(cfg.dataset, 'buffer://odin:1972');
    chanindx = 1:9; % odin buffer specific
else
    chanindx = strmatch('headloc', chantype(hdr));
end

if isempty(chanindx)
    error('the data does not seem to have headlocalization channels');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);

prevSample  = 0;
count       = 0;
UpdatedReference = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true
    
    % determine number of samples available in buffer
    hdr = read_header(cfg.headerfile, 'cache', true);
    
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
        dat = read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here onward it is specific to the headlocalization in the CTF system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % put the data in a fieldtrip-like raw structure
        data.trial{1} = double(dat);
        %data.time{1}  = offset2time(begsample, hdr.Fs, endsample-begsample+1);
        data.time{1}  = (double(begsample) + (0:(double(endsample-begsample))))/hdr.Fs;
        data.label    = hdr.label(chanindx);
        data.hdr      = hdr;
        data.fsample  = hdr.Fs;
                
        x1i = strmatch('HLC0011', data.label);
        y1i = strmatch('HLC0012', data.label);
        z1i = strmatch('HLC0013', data.label);
        x2i = strmatch('HLC0021', data.label);
        y2i = strmatch('HLC0022', data.label);
        z2i = strmatch('HLC0023', data.label);
        x3i = strmatch('HLC0031', data.label);
        y3i = strmatch('HLC0032', data.label);
        z3i = strmatch('HLC0033', data.label);
        
        % convert from meter to cm
        coil1 = data.trial{1}([x1i y1i z1i],:) * 100;
        coil2 = data.trial{1}([x2i y2i z2i],:) * 100;
        coil3 = data.trial{1}([x3i y3i z3i],:) * 100;
        
        figure(1);
        a = subplot(1,2,1);
        h = get(a, 'children');
        hold on;
        
        if ~isempty(h)
            % done on every iteration
            delete(h);
        end
        
        if ~isempty(template)
            % plot the three fiducial positions from the template
            % headcoordinate file
            plot3(template.fid.pnt(1,1), template.fid.pnt(1,2), template.fid.pnt(1,3), 'k^', 'MarkerSize',27,'LineWidth',2);
            plot3(template.fid.pnt(2,1), template.fid.pnt(2,2), template.fid.pnt(2,3), 'ko', 'MarkerSize',27,'LineWidth',2);
            plot3(template.fid.pnt(3,1), template.fid.pnt(3,2), template.fid.pnt(3,3), 'ko', 'MarkerSize',27,'LineWidth',2);
        end
        
        % 'u' keypress updates the reference coordinates
        keypress = get(1,'CurrentCharacter');
        if keypress == 'u'
            fprintf('>>> Updating reference coordinates <<< \n')
            
            % update the reference positions
            UpdatedReference(1,1) = coil1(1,end);
            UpdatedReference(1,2) = coil1(2,end);
            UpdatedReference(1,3) = coil1(3,end);
            UpdatedReference(2,1) = coil2(1,end);
            UpdatedReference(2,2) = coil2(2,end);
            UpdatedReference(2,3) = coil2(3,end);
            UpdatedReference(3,1) = coil3(1,end);
            UpdatedReference(3,2) = coil3(2,end);
            UpdatedReference(3,3) = coil3(3,end);
 
            template = []; % switch off input template
            set(1,'CurrentCharacter','o'); % update switch
        end
        
        if ~isempty(UpdatedReference)
            % plot the updated reference
            plot3(UpdatedReference(1,1),UpdatedReference(1,2),UpdatedReference(1,3), 'k^', 'LineWidth',2,'MarkerSize',27) % nasion
            plot3(UpdatedReference(2,1),UpdatedReference(2,2),UpdatedReference(2,3), 'ko', 'LineWidth',2,'MarkerSize',27) % right ear
            plot3(UpdatedReference(3,1),UpdatedReference(3,2),UpdatedReference(3,3), 'ko', 'LineWidth',2,'MarkerSize',27) % left ear
        end
        
        % plot the coil positons
        plot3(coil1(1,:),coil1(2,:),coil1(3,:), 'k^', 'LineWidth',1,'MarkerSize',3) % nasion
        plot3(coil2(1,:),coil2(2,:),coil2(3,:), 'ko', 'LineWidth',1,'MarkerSize',3) % right ear
        plot3(coil3(1,:),coil3(2,:),coil3(3,:), 'ko', 'LineWidth',1,'MarkerSize',3) % left ear
        
        % compute Center of Mass
        com(1) = (coil2(1,end) + coil3(1,end)) / 2;
        com(2) = (coil2(2,end) + coil3(2,end)) / 2;
        com(3) = (coil2(3,end) + coil3(3,end)) / 2;
        
        % create two vectors from the vertices: nasion, CoM, P (which has nasion_x and _y and CoM_z)
        v1 = [coil1(1,end) - com(1), coil1(2,end) - com(2), coil1(3,end) - com(3)];
        v2 = [coil1(1,end) - com(1), coil1(2,end) - com(2), com(3) - com(3)];
        % find the angle
        theta = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        % convert it to degrees
        angle_nasion_sagittal = (theta * (180/pi));
        
        % check for nasion position
        if ~isempty(UpdatedReference)
            if abs(UpdatedReference(1,1))-cfg.accuracy_green < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(UpdatedReference(1,1))+cfg.accuracy_green ...
                    && abs(UpdatedReference(1,2))-cfg.accuracy_green < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(UpdatedReference(1,2))+cfg.accuracy_green ...
                    && abs(UpdatedReference(1,3))-cfg.accuracy_green < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(UpdatedReference(1,3))+cfg.accuracy_green
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'g^', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
                head1 = true;
            elseif abs(UpdatedReference(1,1))-cfg.accuracy_orange < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(UpdatedReference(1,1))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(1,2))-cfg.accuracy_orange < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(UpdatedReference(1,2))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(1,3))-cfg.accuracy_orange < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(UpdatedReference(1,3))+cfg.accuracy_orange
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'y^', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
                head1 = false;
            else % when not in correct position
                plot3(coil1(1,end),coil1(2,end), coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
                head1 = false;
            end          
        elseif ~isempty(template)
            if abs(template.fid.pnt(1,1))-cfg.accuracy_green < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(template.fid.pnt(1,1))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(1,2))-cfg.accuracy_green < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(template.fid.pnt(1,2))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(1,3))-cfg.accuracy_green < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(template.fid.pnt(1,3))+cfg.accuracy_green
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'g^', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
                head1 = true;
            elseif abs(template.fid.pnt(1,1))-cfg.accuracy_orange < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(template.fid.pnt(1,1))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(1,2))-cfg.accuracy_orange < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(template.fid.pnt(1,2))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(1,3))-cfg.accuracy_orange < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(template.fid.pnt(1,3))+cfg.accuracy_orange
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'y^', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
                head1 = false;
            else % when not in correct position
                plot3(coil1(1,end),coil1(2,end), coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
                head1 = false;
            end
        else
            plot3(coil1(1,end),coil1(2,end), coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            head1 = false;
        end
        
        % check for left ear position
        if ~isempty(UpdatedReference)
            if abs(UpdatedReference(2,1))-cfg.accuracy_green < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(UpdatedReference(2,1))+cfg.accuracy_green ...
                    && abs(UpdatedReference(2,2))-cfg.accuracy_green < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(UpdatedReference(2,2))+cfg.accuracy_green ...
                    && abs(UpdatedReference(2,3))-cfg.accuracy_green < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(UpdatedReference(2,3))+cfg.accuracy_green
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
                head2 = true;
            elseif abs(UpdatedReference(2,1))-cfg.accuracy_orange < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(UpdatedReference(2,1))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(2,2))-cfg.accuracy_orange < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(UpdatedReference(2,2))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(2,3))-cfg.accuracy_orange < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(UpdatedReference(2,3))+cfg.accuracy_orange
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
                head2 = false;
            else % when not in correct position
                plot3(coil2(1,end),coil2(2,end), coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
                head2 = false;
            end
        elseif ~isempty(template)
            if abs(template.fid.pnt(2,1))-cfg.accuracy_green < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(template.fid.pnt(2,1))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(2,2))-cfg.accuracy_green < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(template.fid.pnt(2,2))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(2,3))-cfg.accuracy_green < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(template.fid.pnt(2,3))+cfg.accuracy_green
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
                head2 = true;
            elseif abs(template.fid.pnt(2,1))-cfg.accuracy_orange < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(template.fid.pnt(2,1))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(2,2))-cfg.accuracy_orange < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(template.fid.pnt(2,2))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(2,3))-cfg.accuracy_orange < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(template.fid.pnt(2,3))+cfg.accuracy_orange
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
                head2 = false;
            else % when not in correct position
                plot3(coil2(1,end),coil2(2,end), coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
                head2 = false;
            end
        else
            plot3(coil2(1,end),coil2(2,end), coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            head2 = false;
        end
        
        % check for right ear position
        if ~isempty(UpdatedReference)
            if abs(UpdatedReference(3,1))-cfg.accuracy_green < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(UpdatedReference(3,1))+cfg.accuracy_green  ...
                    && abs(UpdatedReference(3,2))-cfg.accuracy_green  < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(UpdatedReference(3,2))+cfg.accuracy_green  ...
                    && abs(UpdatedReference(3,3))-cfg.accuracy_green  < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(UpdatedReference(3,3))+cfg.accuracy_green
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
                head3 = true;
            elseif abs(UpdatedReference(3,1))-cfg.accuracy_orange < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(UpdatedReference(3,1))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(3,2))-cfg.accuracy_orange < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(UpdatedReference(3,2))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(3,3))-cfg.accuracy_orange < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(UpdatedReference(3,3))+cfg.accuracy_orange
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
                head3 = false;
            else % when not in correct position
                plot3(coil3(1,end),coil3(2,end), coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
                head3 = false;
            end
        elseif ~isempty(template)
            if abs(template.fid.pnt(3,1))-cfg.accuracy_green < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(template.fid.pnt(3,1))+cfg.accuracy_green  ...
                    && abs(template.fid.pnt(3,2))-cfg.accuracy_green  < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(template.fid.pnt(3,2))+cfg.accuracy_green  ...
                    && abs(template.fid.pnt(3,3))-cfg.accuracy_green  < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(template.fid.pnt(3,3))+cfg.accuracy_green
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
                head3 = true;
            elseif abs(template.fid.pnt(3,1))-cfg.accuracy_orange < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(template.fid.pnt(3,1))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(3,2))-cfg.accuracy_orange < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(template.fid.pnt(3,2))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(3,3))-cfg.accuracy_orange < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(template.fid.pnt(3,3))+cfg.accuracy_orange
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
                head3 = false;
            else % when not in correct position
                plot3(coil3(1,end),coil3(2,end), coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
                head3 = false;
            end
        else
            plot3(coil3(1,end),coil3(2,end), coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            head3 = false;
        end
        
        % draw 3d head
        cc = circumcenter(coil1,coil2,coil3);
        x_radius = sqrt((coil2(1,end) - cc(1))^2 + (coil2(2,end) - cc(2))^2);
        y_radius = sqrt((coil3(1,end) - cc(1))^2 + (coil3(2,end) - cc(2))^2);
        [xe, ye, ze] = ellipsoid(cc(1),cc(2),cc(3),x_radius,y_radius,11);
        hh = surfl(xe, ye, ze);
        shading interp
        if head1 == true && head2 == true && head3 == true
            colormap cool
        else
            colormap hot
        end
        alpha(.15)
        
        % plot text
        if ~isempty(template)
            text(-8,8,template.fid.pnt(2,3),'Left','FontSize',15);
            text(6,-6,template.fid.pnt(3,3),'Right','FontSize',15);
        end
        
        % show the update key
        str = sprintf('Press U to update the reference coordinates');
        title(str);
        
        if isempty(h)
            % only done once
            grid on
            if angle_nasion_sagittal > 45
                az = -45;
                el = -90;
            else
                az = -45;
                el = 0;
            end
            xlabel('x (cm)');
            ylabel('y (cm)');
            zlabel('z (cm)');
            set(gca, 'xtick', -30:2.0:30)
            set(gca, 'ytick', -30:2.0:30)
            set(gca, 'ztick', -60:2.0:0) % note the different scaling
            view(az, el)
            % axis square
            axis vis3d
            axis manual
        end
        
        b = subplot(1,2,2);
        i = get(b, 'children');
        hold on;
        
        if ~isempty(i)
            % done on every iteration
            delete(i);
        end
        
        if ~isempty(template)
            % plot the three fiducial positions from the template headcoordinate file
            plot3(template.fid.pnt(1,1), template.fid.pnt(1,2), template.fid.pnt(1,3), 'k^', 'MarkerSize',27,'LineWidth',2);
            plot3(template.fid.pnt(2,1), template.fid.pnt(2,2), template.fid.pnt(2,3), 'ko', 'MarkerSize',27,'LineWidth',2);
            plot3(template.fid.pnt(3,1), template.fid.pnt(3,2), template.fid.pnt(3,3), 'ko', 'MarkerSize',27,'LineWidth',2);
        end
        
        if ~isempty(UpdatedReference)
            % plot the updated reference
            plot3(UpdatedReference(1,1),UpdatedReference(1,2),UpdatedReference(1,3), 'k^', 'LineWidth',2,'MarkerSize',27) % nasion
            plot3(UpdatedReference(2,1),UpdatedReference(2,2),UpdatedReference(2,3), 'ko', 'LineWidth',2,'MarkerSize',27) % right ear
            plot3(UpdatedReference(3,1),UpdatedReference(3,2),UpdatedReference(3,3), 'ko', 'LineWidth',2,'MarkerSize',27) % left ear
        end
        
        % plot the coil positons
        plot3(coil1(1,:),coil1(2,:),coil1(3,:), 'k^', 'LineWidth',1,'MarkerSize',3) % nasion
        plot3(coil2(1,:),coil2(2,:),coil2(3,:), 'ko', 'LineWidth',1,'MarkerSize',3) % right ear
        plot3(coil3(1,:),coil3(2,:),coil3(3,:), 'ko', 'LineWidth',1,'MarkerSize',3) % left ear
        
        % check for nasion position
        if ~isempty(UpdatedReference)
            if abs(UpdatedReference(1,1))-cfg.accuracy_green < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(UpdatedReference(1,1))+cfg.accuracy_green ...
                    && abs(UpdatedReference(1,2))-cfg.accuracy_green < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(UpdatedReference(1,2))+cfg.accuracy_green ...
                    && abs(UpdatedReference(1,3))-cfg.accuracy_green < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(UpdatedReference(1,3))+cfg.accuracy_green
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'g^', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
            elseif abs(UpdatedReference(1,1))-cfg.accuracy_orange < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(UpdatedReference(1,1))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(1,2))-cfg.accuracy_orange < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(UpdatedReference(1,2))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(1,3))-cfg.accuracy_orange < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(UpdatedReference(1,3))+cfg.accuracy_orange
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'y^', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
            else % when not in correct position
                plot3(coil1(1,end),coil1(2,end), coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            end           
        elseif ~isempty(template)
            if abs(template.fid.pnt(1,1))-cfg.accuracy_green < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(template.fid.pnt(1,1))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(1,2))-cfg.accuracy_green < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(template.fid.pnt(1,2))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(1,3))-cfg.accuracy_green < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(template.fid.pnt(1,3))+cfg.accuracy_green
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'g^', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
            elseif abs(template.fid.pnt(1,1))-cfg.accuracy_orange < abs(coil1(1,end)) && abs(coil1(1,end)) < abs(template.fid.pnt(1,1))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(1,2))-cfg.accuracy_orange < abs(coil1(2,end)) && abs(coil1(2,end)) < abs(template.fid.pnt(1,2))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(1,3))-cfg.accuracy_orange < abs(coil1(3,end)) && abs(coil1(3,end)) < abs(template.fid.pnt(1,3))+cfg.accuracy_orange
                plot3(coil1(1,end),coil1(2,end),coil1(3,end),'y^', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
            else % when not in correct position
                plot3(coil1(1,end),coil1(2,end), coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            end
        else
            plot3(coil1(1,end),coil1(2,end), coil1(3,end),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        end
        
        % check for left ear position
        if ~isempty(UpdatedReference)
            if abs(UpdatedReference(2,1))-cfg.accuracy_green < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(UpdatedReference(2,1))+cfg.accuracy_green ...
                    && abs(UpdatedReference(2,2))-cfg.accuracy_green < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(UpdatedReference(2,2))+cfg.accuracy_green ...
                    && abs(UpdatedReference(2,3))-cfg.accuracy_green < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(UpdatedReference(2,3))+cfg.accuracy_green
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
            elseif abs(UpdatedReference(2,1))-cfg.accuracy_orange < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(UpdatedReference(2,1))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(2,2))-cfg.accuracy_orange < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(UpdatedReference(2,2))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(2,3))-cfg.accuracy_orange < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(UpdatedReference(2,3))+cfg.accuracy_orange
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
            else % when not in correct position
                plot3(coil2(1,end),coil2(2,end), coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            end
        elseif ~isempty(template)
            if abs(template.fid.pnt(2,1))-cfg.accuracy_green < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(template.fid.pnt(2,1))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(2,2))-cfg.accuracy_green < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(template.fid.pnt(2,2))+cfg.accuracy_green ...
                    && abs(template.fid.pnt(2,3))-cfg.accuracy_green < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(template.fid.pnt(2,3))+cfg.accuracy_green
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
            elseif abs(template.fid.pnt(2,1))-cfg.accuracy_orange < abs(coil2(1,end)) && abs(coil2(1,end)) < abs(template.fid.pnt(2,1))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(2,2))-cfg.accuracy_orange < abs(coil2(2,end)) && abs(coil2(2,end)) < abs(template.fid.pnt(2,2))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(2,3))-cfg.accuracy_orange < abs(coil2(3,end)) && abs(coil2(3,end)) < abs(template.fid.pnt(2,3))+cfg.accuracy_orange
                plot3(coil2(1,end),coil2(2,end),coil2(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
            else % when not in correct position
                plot3(coil2(1,end),coil2(2,end), coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            end
        else
            plot3(coil2(1,end),coil2(2,end), coil2(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        end
        
        % check for right ear position
        if ~isempty(UpdatedReference)
            if abs(UpdatedReference(3,1))-cfg.accuracy_green < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(UpdatedReference(3,1))+cfg.accuracy_green  ...
                    && abs(UpdatedReference(3,2))-cfg.accuracy_green  < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(UpdatedReference(3,2))+cfg.accuracy_green  ...
                    && abs(UpdatedReference(3,3))-cfg.accuracy_green  < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(UpdatedReference(3,3))+cfg.accuracy_green
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
            elseif abs(UpdatedReference(3,1))-cfg.accuracy_orange < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(UpdatedReference(3,1))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(3,2))-cfg.accuracy_orange < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(UpdatedReference(3,2))+cfg.accuracy_orange ...
                    && abs(UpdatedReference(3,3))-cfg.accuracy_orange < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(UpdatedReference(3,3))+cfg.accuracy_orange
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
            else % when not in correct position
                plot3(coil3(1,end),coil3(2,end), coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            end
        elseif ~isempty(template)
            if abs(template.fid.pnt(3,1))-cfg.accuracy_green < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(template.fid.pnt(3,1))+cfg.accuracy_green  ...
                    && abs(template.fid.pnt(3,2))-cfg.accuracy_green  < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(template.fid.pnt(3,2))+cfg.accuracy_green  ...
                    && abs(template.fid.pnt(3,3))-cfg.accuracy_green  < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(template.fid.pnt(3,3))+cfg.accuracy_green
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
            elseif abs(template.fid.pnt(3,1))-cfg.accuracy_orange < abs(coil3(1,end)) && abs(coil3(1,end)) < abs(template.fid.pnt(3,1))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(3,2))-cfg.accuracy_orange < abs(coil3(2,end)) && abs(coil3(2,end)) < abs(template.fid.pnt(3,2))+cfg.accuracy_orange ...
                    && abs(template.fid.pnt(3,3))-cfg.accuracy_orange < abs(coil3(3,end)) && abs(coil3(3,end)) < abs(template.fid.pnt(3,3))+cfg.accuracy_orange
                plot3(coil3(1,end),coil3(2,end),coil3(3,end),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
            else % when not in correct position
                plot3(coil3(1,end),coil3(2,end), coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
            end
        else
            plot3(coil3(1,end),coil3(2,end), coil3(3,end),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        end
        
        % plot head
        ih = surfl(xe, ye, ze);
        shading interp
        if head1 == true && head2 == true && head3 == true
            colormap cool
        else
            colormap hot
        end
        alpha(.15)
        
        % plot text
        if ~isempty(template)
            text(-8,8,template.fid.pnt(2,3),'Left','FontSize',15);
            text(6,-6,template.fid.pnt(3,3),'Right','FontSize',15);
        end
        
        % show current timesample
        str = sprintf('Press CTRL+C to stop, time = %d s\n', round(mean(data.time{1})));
        title(str);
        fprintf(str);
        
        if isempty(i)
            % only done once
            grid on
            if angle_nasion_sagittal > 45
                az = -45;
                el = 0;
            else
                az = -45;
                el = 90;
            end
            xlabel('x (cm)');
            ylabel('y (cm)');
            zlabel('z (cm)');
            set(gca, 'xtick', -30:2.0:30)
            set(gca, 'ytick', -30:2.0:30)
            set(gca, 'ztick', -60:2.0:0) % note the different scaling
            view(az, el)
            % axis square
            axis vis3d
            axis square
        end
        
        % force Matlab to update the figure
        drawnow
        
    end % if enough new samples
end % while true


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cc] = circumcenter(coil1,coil2,coil3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the circumcenter of our 3D triangle (3 coils)
% output: cc(x,y,z)

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
