function nmt_spmfig_setup(cfg)
% adds custom features to SPM MRI display window
% e.g., head coords, MNI coords, activation intensity
% Author: Sarang S. Dalal, NEMOlab

global st

fg = spm_figure('GetWin','Graphics');
WS = spm('winscale');

% reposition MRI intensity box
inlabel = findobj('String','Intensity:','HorizontalAlignment','center');

% get rid of some stuff, sweep other stuff under the rug
delete(findobj('String','Crosshair Position'));
delete(findobj('ToolTipString','move crosshairs to origin'));
delete(inlabel);
set(st.in,'Visible','off');

beaminlabel = uicontrol(fg,'Style','Text','Position',[60 350 120 020].*WS,'String','Functional:');
st.nmt.beamin = uicontrol(fg,'Style','edit', 'Position',[175 350  125 020].*WS,'String','');

uicontrol(fg,'Style','Text', 'Position',[75 255 35 020].*WS,'String','MNI:');
% st.mnip = uicontrol(fg,'Style','edit', 'Position',[110 255 135 020].*WS,'String','','Callback','nmt_image(''setposmni'')','ToolTipString','move crosshairs to MNI mm coordinates');
% st.mnilabel=uicontrol('Style','text','BackgroundColor',[1 1 1],'Units','normalized','Position',[.6 .5 .3 .15],'FontSize',12);

uicontrol(fg,'Style','Text', 'Position',[75 315 35 020].*WS,'String','MEG:');
% st.megp = uicontrol(fg,'Style','edit', 'Position',[110 315 135 020].*WS,'String',sprintf('%.1f %.1f %.1f',nut_mri2meg(spm_orthviews('pos')')),'Callback','nmt_image(''setposmeg'')','ToolTipString','move crosshairs to MEG mm coordinates');
st.nmt.megp = uicontrol(fg,'Style','edit', 'Position',[110 315 135 020].*WS,'String',sprintf('%.1f %.1f %.1f',(spm_orthviews('pos')')),'Callback','','ToolTipString','move crosshairs to MEG mm coordinates');

set(st.mp,'Callback','spm_image(''setposmm''); nmt_image(''shopos'');');
set(st.vp,'Callback','spm_image(''setposvx''); nmt_image(''shopos'');');

for ii=1:7
   spmUIhL = findobj('Callback',['spm_image(''repos'',' num2str(ii) ')']);
   set(spmUIhL,'Callback',['spm_image(''repos'',' num2str(ii) '); nmt_image']);
end

spmUIhR = [findobj('Style','popupmenu'); findobj('Style','togglemenu')];
for ii=1:length(spmUIhR)
    cb = get(spmUIhR(ii),'Callback');
    set(spmUIhR(ii),'Callback',[cb '; nmt_image;']);
end

% %% so instead, just destroy them all
% hitlist = ...
% [75 220 100 016
% 75 200 100 016
% 75 180 100 016
% 75 160 100 016 
% 75 140 100 016 
% 75 120 100 016 
% 75 100 100 016 
% 75  80 100 016 
% 75  60 100 016 
% 175 220 065 020 
% 175 200 065 020 
% 175 180 065 020 
% 175 160 065 020 
% 175 140 065 020 
% 175 120 065 020 
% 175 100 065 020 
% 175  80 065 020 
% 175  60 065 020
% 70 35 125 020
% 195 35 55 020
% ];
% 
% for ii=1:size(hitlist,1)
%     delete(findobj('Position',hitlist(ii,:).*WS));
% end
%     




%% add sliders for scrolling through MRI
% if ndefaults.sliders
%     ax_pos = get(st.vols{1}.ax{1}.ax,'Position');
%     cor_pos = get(st.vols{1}.ax{2}.ax,'Position');
%     sag_pos = get(st.vols{1}.ax{3}.ax,'Position');
%     slidercallback1 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([st.vols{1}.dim(1)*get(gcbo,''Value'') pos(2) pos(3)])); spm_image(''setposvx''); nmt_image(''shopos'')';
%     slidercallback2 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([pos(1) st.vols{1}.dim(2)*get(gcbo,''Value'') pos(3)])); spm_image(''setposvx''); nmt_image(''shopos'')';
%     slidercallback3 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([pos(1) pos(2) st.vols{1}.dim(3)*get(gcbo,''Value'')])); spm_image(''setposvx''); nmt_image(''shopos'')';
%     
%     mat = st.vols{1}.premul*st.vols{1}.mat;
%     R=spm_imatrix(st.vols{1}.mat);
%     R = spm_matrix([0 0 0 R(4:6)]);
%     R = R(1:3,1:3);
%     dim = st.vols{1}.dim*inv(R);
%     dim_mm = abs(dim .* diag(st.vols{1}.mat(1:3,1:3)*inv(R))');
%     
%     st.slider{1}=uicontrol(fg,'Style','Slider','Callback',slidercallback1,'SliderStep',[voxelsize(1)/abs(dim_mm(1)) .1],'Units','normalized','Position',[ax_pos(1) ax_pos(2)+ax_pos(4) ax_pos(3) .015]);
%     st.slider{2}=uicontrol(fg,'Style','Slider','Callback',slidercallback2,'SliderStep',[voxelsize(2)/dim_mm(2) .1],'Units','normalized','Position',[ax_pos(1)+ax_pos(3) ax_pos(2) .02 ax_pos(4)]);
%     st.slider{3}=uicontrol(fg,'Style','Slider','Callback',slidercallback3,'SliderStep',[voxelsize(3)/dim_mm(3) .1],'Units','normalized','Position',[cor_pos(1)+cor_pos(3) cor_pos(2) .02 cor_pos(4)]);
%     set([st.slider{:}],'Visible',rivets.sliderenable);
% end


for i=1:3  % step through three orthogonal views
    SPM_axes_obj(i) = st.vols{1}.ax{i}.ax;
end

% for time series or spectrogram data, expand SPM window and create new axes
if(isfield(st.nmt,''))
end

switch(cfg.funparameter)
    case {'avg.mom','mom'}
        % MRI slices have "relative" position; change to fixed position
        un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');

        set(SPM_axes_obj,'Units','pixels');
        % now we can expand SPM8 window, and MRI stays in place
        winpos = get(fg,'Position');
        set(fg,'Position',[1 40 2*winpos(3) winpos(4)]);
        
        set(SPM_axes_obj,'Units','normalized');
        
        for ii=1:3
            st.vols{1}.ax{ii}.axpos = get(SPM_axes_obj(ii),'Position');
        end
        
        %% modeled after section in spm_orthviews.m that sets axis positions
        Dims = diff(st.bb)'+1;
        
        sz    = get(st.fig,'Position');set(st.fig,'Units',un);
        sz    = sz(3:4);
        sz(2) = sz(2)-40;
        
        area = st.vols{1}.area(:);
        
        area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
        if st.mode == 0,
            sx   = area(3)/(Dims(1)+Dims(3))/1.02;
        else
            sx   = area(3)/(Dims(1)+Dims(2))/1.02;
        end;
        sy   = area(4)/(Dims(2)+Dims(3))/1.02;
        s    = min([sx sy]);
        
        offy = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
        sky = s*(Dims(2)+Dims(3))*0.02;
        if st.mode == 0,
            offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
            skx = s*(Dims(1)+Dims(3))*0.02;
        else
            offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
            skx = s*(Dims(1)+Dims(2))*0.02;
        end;
        
        %%
        st.nmt.ax_ts = axes('Visible','off','DrawMode','fast','Parent',st.fig,...
            'YDir','normal');
        
        % Transverse
        set(st.nmt.ax_ts,'Units','pixels', ...
            'Position',[offx+s*Dims(1)+4*skx offy+s*Dims(2) s*(Dims(1)+Dims(2)) s*(Dims(2))],...
            'Units','normalized','Visible','on')
        
        % save new MRI axis positions, because original SPM controls may mess them up
        for ii=1:3
            st.nmt.mriaxpos(:,ii) = get(st.vols{1}.ax{ii}.ax,'Position');
        end
end

tmp = get(SPM_axes_obj(1),'ButtonDownFcn');
if(isa(tmp,'function_handle')) % this happens only with raw SPM8 ButtonDownFcn
    st.nmt.spm_refresh_handle = tmp;
end
set(SPM_axes_obj,'ButtonDownFcn',('nmt_image(''shopos'')'));

[~,mrifile]=fileparts(cfg.mripath);
% if( isfield(coreg,'norm_mripath') || strcmp(mrifile(1),'w') || strcmp(mrifile,'T1') || strncmp(mrifile,'avg',3) )
%     rivets.TalDB = load('talairachDB');
% end

