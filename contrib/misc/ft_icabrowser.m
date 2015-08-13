function [rej_comp] = ft_icabrowser(cfg, comp)

% ICA component viewer and GUI
%
% loads in comp structure from FieldTrip ft_componentanalysis
% presents a GUI interface showin the power spectrum, variance over time
% and the topography of the components, as well as the possibility to save
% a PDF, view the timecourse and toggle components to be rejected vs kept.
% when done, will create a file with the components to be rejected
%
% CONFIGURATION NEEDED:
% cfg.path         where pdfs will be saves
% cfg.prefix       prefix of the pdf files
% cfg.layout       layout of the topo view
%
% OPTIONAL CONFIGURATION:
% cfg.colormap      colormap for topo
% cfg.inputfile
% cfg.outputfile    will contain indices of all components to reject
%
% original written by Thomas Pfeffer
% adapted by Jonathan Daume and Anne Urai
% University Medical Center Hamburg-Eppendorf, 2015

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(cfg, 'inputfile'), load(cfg.inputfile); end
assert(exist('comp', 'var') > 0, 'could not find a comp structure');

% setup
var_data = cat(2,comp.trial{:});
var_time = (1:(size(var_data,2)))/comp.fsample;
% only do the fft on a subset of trials, saves time
fft_data = cat(2,comp.trial{1:5:end});
% preallocate rejected components
rej_comp = zeros(size(comp.label,1),1);

subpl = 4;
l = 0;
cnt = 1;
il = 0;
logar = 1; % default logarithmic or linear scale for powspctrm

% set path
path = cfg.path;
if ~exist(path, 'dir'),
    mkdir path;
end
prefix = cfg.prefix;

% to save time redoing this for each topo
cfglay = keepfields(cfg, {'layout'});
lay = ft_prepare_layout(cfglay, comp);

cfgtopo = [];
cfgtopo.layout    = lay;     % specify the layout file that should be used for plotting
cfgtopo.comment   = 'no';
cfgtopo.highlight = 'off';
cfgtopo.marker    = 'off';
cfgtopo.style     = 'straight';
if isfield(cfg, 'colormap'),  cfgtopo.colormap  = cfg.colormap;  end

err = 0;
manpos = [0.1 0.1 0.8 0.8]; % figure position, can be updated later

% ------------------------------------------------
% COMPUTE LATENCY FOR 2s-WINDOWS
% ------------------------------------------------

slen = floor(2*comp.fsample);
smax = floor(size(var_data,2)/slen);
comp_var  = nan(subpl, smax); % preallocate
comp_time = nan(1, smax);     % preallocate
for s = 1 : smax
  comp_time(s) = mean(var_time(1,(s-1)*slen+1:s*slen));
end
        
while err == 0 % KEEP GOING UNTIL THERE IS AN ERROR
    
    while il < subpl, % il is the subplot count
        
        il = il + 1;
        i = (cnt-1)*subpl+il;
        
        if mod(i-1,subpl)==0
            % keep manual screen position - better in dual monitor settings
            f = figure('units','normalized','outerposition', manpos);
            l = l + 1;
        end
        
        % ------------------------------------------------
        % COMPUTE VARIANCE FOR 2s-WINDOWS
        % ------------------------------------------------
        
        for s = 1 : smax
            comp_var(i,s)=var(var_data(i,(s-1)*slen+1:s*slen));
        end
        
        % ------------------------------------------------
        % COMPUTE POWER SPECTRUM
        % ------------------------------------------------
        smo = 50;
        steps = 10;
        Fs = comp.fsample;
        N = floor(size(fft_data,2));
        xdft = fft(fft_data(i,:));
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)).*abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        
        j = 1;
        k = 1;
        while j < length(psdx)-smo
            smoothed(k)=mean(psdx(j:j+smo));
            j = j + steps;
            k = k + 1;
        end
        
        freq = linspace(0,Fs/2,size(smoothed,2));
        strt = find(freq > 2,1,'first');
        stp  = find(freq < 200,1,'last');
        
        % ------------------------------------------------
        % PLOT POWER SPECTRUM
        % ------------------------------------------------
        subcomp{1}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-2);
        if logar
            plot(freq(strt:stp),log10(smoothed(strt:stp)));
            ylabel('(dB/Hz)');
        else
            plot(freq(strt:stp),smoothed(strt:stp));
            ylabel('T^2/Hz');
        end
        set(gca,'TickDir','out','XTick',0:25:200)
        xlabel('Frequency (Hz)'); grid on;
        axis tight;
        
        % ------------------------------------------------
        % PLOT VARIANCE OVER TIME
        % ------------------------------------------------
        subcomp{2}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-1);
        scatter(comp_time,comp_var(i,:),'k.');
        xlabel('Time (s)'); ylabel('Variance');
        axis tight; set(gca, 'tickdir', 'out');
        
        % ------------------------------------------------
        % PLOT COMPONENT TOPOGRAPHY
        % ------------------------------------------------
        subcomp{3}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3);
        cfgtopo.component = i;       % specify the component(s) that should be plotted
        ft_topoplotIC(cfgtopo, comp);
        
        if mod(i,subpl)==0 || i == 80
            
            pos = [0.76 0.73 0.075 0.035; ...
                0.76 0.51 0.075 0.035; ...
                0.76 0.29 0.075 0.035; ...
                0.76 0.07 0.075 0.035];
            
            for ibgc = 1 : 4
                if rej_comp(i+(ibgc-4)) == 1
                    bgc{ibgc} = 'r';
                else
                    bgc{ibgc} = 'g';
                end
            end
            
            % ------------------------------------------------
            % SHOW TIMECOURSE OF THIS COMPONENT
            % ------------------------------------------------
            tc1 = uicontrol('Units','normalized','Position',[0.86 0.73 0.075 0.035],...
                'Style','pushbutton','String','Timecourse','Callback',{@tc_cb, 1});
            tc2 = uicontrol('Units','normalized','Position',[0.86 0.51 0.075 0.035],...
                'Style','pushbutton','String','Timecourse','Callback',{@tc_cb, 2});
            tc3 = uicontrol('Units','normalized','Position',[0.86 0.29 0.075 0.035],...
                'Style','pushbutton','String','Timecourse','Callback',{@tc_cb, 3});
            tc4 = uicontrol('Units','normalized','Position',[0.86 0.07 0.075 0.035],...
                'Style','pushbutton','String','Timecourse','Callback',{@tc_cb, 4});
            
            % ------------------------------------------------
            % REJECT COMPONENT
            % ------------------------------------------------
            rej1 = uicontrol('Units','normalized', 'Tag', 'rej1', 'Position',pos(1,:),'Style','pushbutton','String','Keep', ...
                'Backgroundcolor',bgc{1},'Callback',@rej_callback1);
            rej2 = uicontrol('Units','normalized','Tag', 'rej2','Position',pos(2,:),'Style','pushbutton','String','Keep', ...
                'Backgroundcolor',bgc{2},'Callback',@rej_callback2);
            rej3 = uicontrol('Units','normalized','Tag', 'rej3','Position',pos(3,:),'Style','pushbutton','String','Keep', ...
                'Backgroundcolor',bgc{3},'Callback',@rej_callback3);
            rej4 = uicontrol('Units','normalized','Tag', 'rej4','Position',pos(4,:),'Style','pushbutton','String','Keep', ...
                'Backgroundcolor',bgc{4},'Callback',@rej_callback4);
            
            % ------------------------------------------------
            % SAVE COMPONENT PDF
            % ------------------------------------------------
            savecomp1 = uicontrol('Units','normalized','Position',[0.86 0.78 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',{@sc_cb, 1});
            savecomp2 = uicontrol('Units','normalized','Position',[0.86 0.56 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',{@sc_cb, 2});
            savecomp3 = uicontrol('Units','normalized','Position',[0.86 0.34 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',{@sc_cb, 3});
            savecomp4 = uicontrol('Units','normalized','Position',[0.86 0.12 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',{@sc_cb, 4});
            
            % ------------------------------------------------
            % MOVE TO NEXT
            % ------------------------------------------------
            
            if i > 4
                prev = uicontrol('Units','normalized','Position',[0.1 0.01 0.075 0.05],'Style','pushbutton','String','Prev','Callback',@nextplot);
            else
                prev = uicontrol('Units','normalized','Position',[0.1 0.01 0.075 0.05],'Style','pushbutton','String','');
            end
            
            if i < size(comp.label,1)-3
                next = uicontrol('Units','normalized','Position',[0.2 0.01 0.075 0.05],'Style','pushbutton','String','Next','Callback',@lastplot);
            else
                next = uicontrol('Units','normalized','Position',[0.2 0.01 0.075 0.05],'Style','pushbutton','String','');
            end
            
            % ------------------------------------------------
            % CHANGE FREQUENCY AXIS
            % ------------------------------------------------
            logarith = uicontrol('Units','normalized','Position',[0.3 0.01 0.075 0.05],'Style','pushbutton','String','Log10','Callback',@plotlog);
            lin      = uicontrol('Units','normalized','Position',[0.4 0.01 0.075 0.05],'Style','pushbutton','String','Linear','Callback',@plotlin);
            
            % ------------------------------------------------
            % SAVE AND QUIT
            % ------------------------------------------------
            s = uicontrol('Units','normalized','Position',[0.90 0.01 0.075 0.05], ...
                'Style','pushbutton','String','Save','Callback',@save_callback);
            quit_it = uicontrol('Units','normalized','Position',[0.80 0.01 0.075 0.05],...
                'Style','pushbutton','String','Quit','Callback',@quitme);
            orient landscape
            uiwait
        end
    end
end

% ------------------------------------------------
% DEFINE NESTED CALLBACK FUNCTIONS
% ------------------------------------------------
    function rej_callback1(h, evt)
        rej1 = findobj('Tag','rej1');
        if (rej_comp(i-3) == 0),
            set(rej1,'Backgroundcolor','r'),
            rej_comp(i-3)=1;
        else
            set(rej1,'Backgroundcolor','g'),
            rej_comp(i-3)=0;
        end
    end

    function rej_callback2(h, evt)
        rej2 = findobj('Tag','rej2');
        if (rej_comp(i-2) == 0),
            set(rej2,'Backgroundcolor','r'),
            rej_comp(i-2)=1;
        else
            set(rej2,'Backgroundcolor','g'),
            rej_comp(i-2)=0;
        end
    end
    function rej_callback3(h, evt)
        rej3 = findobj('Tag','rej3');
        if (rej_comp(i-1) == 0),
            set(rej3,'Backgroundcolor','r'),
            rej_comp(i-1)=1;
        else
            set(rej3,'Backgroundcolor','g'),
            rej_comp(i-1)=0;
        end
    end
    function rej_callback4(h, evt)
        rej4 = findobj('Tag','rej4');
        if (rej_comp(i-0) == 0),
            set(rej4,'Backgroundcolor','r'),
            rej_comp(i-0)=1;
        else
            set(rej4,'Backgroundcolor','g'),
            rej_comp(i-0)=0;
        end
    end

% timecourse funcs
    function tc_cb(h, evt, whichcomp)
        cfgtc = [];
        cfgtc.layout = lay;
        cfgtc.viewmode = 'butterfly';
        cfgtc.channel = [i+(4-whichcomp)];
        % cfgtc.ylim = [-5e-13 5e-13];
        ft_databrowser(cfgtc, comp);
    end

% save to figure
    function sc_cb(h, evt, whichcomp)
        h = figure;
        set(h,'Position',[200 200 1000 300]);
        set(h,'Units','inches');
        screenposition = get(h,'Position');
        set(h, 'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
        new = copyobj(subcomp{1}{whichcomp},h);
        set(new,'Position',[.05 .1 0.25 0.85]);
        new = copyobj(subcomp{2}{whichcomp},h);
        set(new,'Position',[.35 .1 0.25 0.85]);
        set(new,'LineWidth',2);
        new = copyobj(subcomp{3}{whichcomp},h); set(new,'Position',[.55 .05 0.5 0.95]);
        
        % save under the correct comp nr
        compnrs = [i-3 : i];
        print(h,'-dpdf',sprintf('%s/%s_comp%d.pdf', path, prefix, compnrs(whichcomp)));
        fprintf('saved pdf to %s/%s_comp%d.pdf', path, prefix, compnrs(whichcomp));
        close(h)
    end

% gui
    function nextplot(h, evt)
        manpos = get(f,'OuterPosition');
        cnt = cnt - 1;
        il = 0;
        l = l - 2;
        close;
    end

    function lastplot(h, evt)
        manpos = get(f,'OuterPosition');
        cnt = cnt + 1;
        il = 0; close;
    end

    function plotlog(h, evt)
        manpos = get(f,'OuterPosition');
        il = 0;
        logar = 1;
        l = l - 1; close;
    end

    function plotlin(h, evt)
        manpos = get(f,'OuterPosition');
        il = 0; logar = 0; l = l - 1; close;
    end

    function save_callback(h, evt)
        idx = find(rej_comp==1);
        if isfield(cfg, 'outputfile'),
            save(cfg.outputfile,'idx','rej_comp');
        else
            save(sprintf('%s/%s_rejectedcomps.mat', cfg.path, cfg.prefix), 'idx', 'rej_comp');
        end
        close(f);
        disp('saved');
        err = 1;
    end

    function quitme(h, evt)
        close(f);
        err = 1;
    end

end % main function
