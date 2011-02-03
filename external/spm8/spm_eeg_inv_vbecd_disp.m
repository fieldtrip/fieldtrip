function spm_eeg_inv_vbecd_disp(action,varargin)
% Display the dipoles as obtained from VB-ECD
%
% FORMAT spm_eeg_inv_vbecd_disp('Init',D)
% Display the latest VB-ECD solution saved in the .inv{} field of the
% data structure D.
%
% FORMAT spm_eeg_inv_vbecd_disp('Init',D, ind)
% Display the ind^th .inv{} cell element, if it is actually a VB-ECD 
% solution.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id$

% Note:
% unfortunately I cannot see how to ensure that when zooming in the image
% the dipole location stays in place...

global st

Fig     = spm_figure('GetWin','Graphics');
colors  = {'y','b','g','r','c','m'};              % 6 possible colors
marker  = {'o','x','+','*','s','d','v','p','h'};  % 9 possible markers
Ncolors = length(colors);
Nmarker = length(marker);

if nargin == 0, action = 'Init'; end;

switch lower(action),
    
%==========================================================================
case 'init'
%==========================================================================
% FORMAT spm_eeg_inv_vbecd_disp('init',D,ind)
% Initialise the variables with GUI
%--------------------------------------------------------------------------

if nargin<2
    D = spm_eeg_load;
else
    D = varargin{1};
end
if nargin<3
    % find the latest inverse produced with vbecd
    Ninv = length(D.inv);
    lind = [];
    for ii=1:Ninv
        if isfield(D.inv{ii},'method') && ...
                strcmp(D.inv{ii}.method,'vbecd')
            lind = [lind ii];
        end
    end
    ind = max(lind);
    if ~ind, 
        spm('alert*','No VB-ECD solution found with this data file!',...
                'VB-ECD display')
        return
    end
else
    ind = varargin{3};
end


% Stash dipole(s) information in sdip structure
sdip = D.inv{ind}.inverse;

% if the exit flag is not in the structure, assume everything went ok.
if ~isfield(sdip,'exitflag')
    sdip.exitflag = ones(1,sdip.n_seeds);
end

try
    Pimg = spm_vol(D.inv{ind}.mesh.sMRI);
catch
    Pimg = spm_vol(fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii'));
end

spm_orthviews('Reset');
spm_orthviews('Image', Pimg, [0.0 0.45 1 0.55]);
spm_orthviews('MaxBB');
spm_orthviews('AddContext')
st.callback = 'spm_image(''shopos'');';
% remove clicking in image
for ii=1:3,
    set(st.vols{1}.ax{ii}.ax,'ButtonDownFcn',';');
end
WS = spm('WinScale');

% Build GUI
%==========================================================================
% Location:
%--------------------------------------------------------------------------
uicontrol(Fig,'Style','Frame','Position',[60 25 200 325].*WS, ...
    'DeleteFcn','spm_image(''reset'');');
uicontrol(Fig,'Style','Frame','Position',[70 250 180 90].*WS);
uicontrol(Fig,'Style','Text', 'Position',[75 320 170 016].*WS, ...
    'String','Current Position');
uicontrol(Fig,'Style','Text', 'Position',[75 295 35 020].*WS,'String','mm:');
uicontrol(Fig,'Style','Text', 'Position',[75 275 35 020].*WS,'String','vx:');
uicontrol(Fig,'Style','Text', 'Position',[75 255 75 020].*WS,'String','Img Intens.:');

st.mp = uicontrol(Fig,'Style','Text', 'Position',[110 295 135 020].*WS,'String','');
st.vp = uicontrol(Fig,'Style','Text', 'Position',[110 275 135 020].*WS,'String','');
st.in = uicontrol(Fig,'Style','Text', 'Position',[150 255  85 020].*WS,'String','');

c = 'if get(gco,''Value'')==1, spm_orthviews(''Xhairs'',''off''), else, spm_orthviews(''Xhairs'',''on''); end;';
uicontrol(Fig,'Style','togglebutton','Position',[95 220 125 20].*WS,...
    'String','Hide Crosshairs','Callback',c,'ToolTipString','show/hide crosshairs');

% Dipoles/seeds selection:
%--------------------------------------------------------------------------
uicontrol(Fig,'Style','Frame','Position',[300 25 180 325].*WS);
sdip.hdl.hcl = uicontrol(Fig,'Style','pushbutton','Position',[310 320 100 20].*WS, ...
        'String','Clear all','CallBack','spm_eeg_inv_vbecd_disp(''ClearAll'')');

sdip.hdl.hseed=zeros(sdip.n_seeds,1);
for ii=1:sdip.n_seeds
    if sdip.exitflag(ii)==1
        sdip.hdl.hseed(ii) = uicontrol(Fig,'Style','togglebutton','String',num2str(ii),...
            'Position',[310+rem(ii-1,8)*20 295-fix((ii-1)/8)*20 20 20].*WS,...
            'CallBack','spm_eeg_inv_vbecd_disp(''ChgSeed'')');
    else
        sdip.hdl.hseed(ii) = uicontrol(Fig,'Style','Text','String',num2str(ii), ...
            'Position',[310+rem(ii-1,8)*20 293-fix((ii-1)/8)*20 20 20].*WS) ;
    end
end

uicontrol(Fig,'Style','text','String','Select dipole # :', ...
        'Position',[310 255-fix((sdip.n_seeds-1)/8)*20 110 20].*WS);

txt_box = cell(sdip.n_dip,1);
for ii=1:sdip.n_dip, txt_box{ii} = num2str(ii); end
txt_box{sdip.n_dip+1} = 'all';
sdip.hdl.hdip = uicontrol(Fig,'Style','popup','String',txt_box, ...
        'Position',[420 258-fix((sdip.n_seeds-1)/8)*20 40 20].*WS, ...
        'Callback','spm_eeg_inv_vbecd_disp(''ChgDip'')');

% Dipoles orientation and strength:
%--------------------------------------------------------------------------
uicontrol(Fig,'Style','Frame','Position',[70 120 180 90].*WS);
uicontrol(Fig,'Style','Text', 'Position',[75 190 170 016].*WS, ...
                'String','Dipole orientation & strength');
uicontrol(Fig,'Style','Text', 'Position',[75 165 65 020].*WS, ...
                'String','x-y-z or.:');
uicontrol(Fig,'Style','Text', 'Position',[75 145 75 020].*WS, ...
                'String','theta-phi or.:');
uicontrol(Fig,'Style','Text', 'Position',[75 125 75 020].*WS, ...
                'String','Dip. intens.:');
   
sdip.hdl.hor1 = uicontrol(Fig,'Style','Text', 'Position', ...
                            [140 165  105 020].*WS,'String','a');
sdip.hdl.hor2 = uicontrol(Fig,'Style','Text', 'Position', ...
                            [150 145  85 020].*WS,'String','b');
sdip.hdl.int  = uicontrol(Fig,'Style','Text', 'Position', ...
                            [150 125  85 020].*WS,'String','c');
    
st.vols{1}.sdip = sdip;

% First plot = all the seeds that converged !
l_conv = find(sdip.exitflag==1);
if isempty(l_conv)
    error('No seed converged towards a stable solution, nothing to be displayed !')
else
    spm_eeg_inv_vbecd_disp('DrawDip',l_conv,1)
    set(sdip.hdl.hseed(l_conv),'Value',1); % toggle all buttons
end


%==========================================================================
case 'drawdip'
%==========================================================================
% FORMAT spm_eeg_inv_vbecd_disp('DrawDip',i_seed,i_dip,sdip)
% e.g. spm_eeg_inv_vbecd_disp('DrawDip',1,1,sdip)
% e.g. spm_eeg_inv_vbecd_disp('DrawDip',[1:5],1,sdip)
%--------------------------------------------------------------------------

if nargin < 2
    i_seed = 1;
else
    i_seed = varargin{1};
end
if nargin<3
    i_dip = 1;
else
    i_dip = varargin{2}; 
end

if nargin<4
    if isfield(st.vols{1},'sdip')
        sdip = st.vols{1}.sdip;
    else
        error('I can''t find sdip structure');
    end
else
    sdip = varargin{3};
    st.vols{1}.sdip = sdip;
end

if any(i_seed>sdip.n_seeds) || i_dip>(sdip.n_dip+1)
    error('Wrong i_seed or i_dip index in spm_eeg_inv_vbecd_disp');
end

% Note if i_dip==(sdip.n_dip+1) all dipoles are displayed simultaneously,
% The 3D cut will then be at the mean location of all sources !!!
if i_dip == (sdip.n_dip+1)
    i_dip = 1:sdip.n_dip;
end

% if seed indexes passed is wrong (no convergence) remove the wrong ones
i_seed(sdip.exitflag(i_seed)~=1) = [];
if isempty(i_seed)
    error('You passed the wrong seed indexes...')
end
if size(i_seed,2)==1, i_seed=i_seed'; end

% Display business
%--------------------------------------------------------------------------
loc_mm = sdip.loc{i_seed(1)}(:,i_dip);
if length(i_seed)>1
%     unit = ones(1,sdip.n_dip);
    for ii = i_seed(2:end)
        loc_mm = loc_mm + sdip.loc{ii}(:,i_dip);
    end
    loc_mm = loc_mm/length(i_seed);
end
if length(i_dip)>1
    loc_mm = mean(loc_mm,2);
end

% Place the underlying image at right cuts
spm_orthviews('Reposition',loc_mm);

if length(i_dip)>1
    tabl_seed_dip = [kron(ones(length(i_dip),1),i_seed') ...
                        kron(i_dip',ones(length(i_seed),1))];
else
    tabl_seed_dip = [i_seed' ones(length(i_seed),1)*i_dip];
end

% Scaling, according to all dipoles in the selected seed sets.
% The time displayed is the one corresponding to the maximum EEG power !
Mn_j = -1;
l3 = -2:0;
for ii = 1:length(i_seed)
    for jj = 1:sdip.n_dip
        Mn_j = max([Mn_j sqrt(sum(sdip.j{ii}(jj*3+l3,sdip.Mtb).^2))]);
    end
end
st.vols{1}.sdip.tabl_seed_dip = tabl_seed_dip;

% Display all dipoles, the 1st one + the ones close enough.
% Run through the 6 colors and 9 markers to differentiate the dipoles.
% NOTA: 2 dipoles coming from the same set will have same colour/marker
ind = 1 ;
dip_h = zeros(9,size(tabl_seed_dip,1),1);
    % each dipole displayed has 9 handles:
    %   3 per view (2*3): for the line, for the circle & for the error
js_m = zeros(3,1);

% Deal with case of multiple i_seed and i_dip displayed.
% make sure dipole from same i_seed have same colour but different marker.

pi_dip = find(diff(tabl_seed_dip(:,2)));
if isempty(pi_dip)
    % i.e. only one dip displayed per seed, use old fashion
    for ii=1:size(tabl_seed_dip,1)
        if ii>1
            if tabl_seed_dip(ii,1)~=tabl_seed_dip(ii-1,1)
                ind = ind+1;
            end
        end
        ic = mod(ind-1,Ncolors)+1;
        im = fix(ind/Ncolors)+1;

        loc_pl = sdip.loc{tabl_seed_dip(ii,1)}(:,tabl_seed_dip(ii,2));
        js = sdip.j{tabl_seed_dip(ii,1)}(tabl_seed_dip(ii,2)*3+l3,sdip.Mtb);
        vloc = sdip.cov_loc{tabl_seed_dip(ii,1)}(tabl_seed_dip(ii,2)*3+l3,tabl_seed_dip(ii,2)*3+l3);
        dip_h(:,ii) = add1dip(loc_pl,js/Mn_j*20,vloc, ...
                            marker{im},colors{ic},st.vols{1}.ax,Fig,st.bb);
        js_m = js_m+js;
    end
else
    for ii=1:pi_dip(1)
        if ii>1
            if tabl_seed_dip(ii,1)~=tabl_seed_dip(ii-1,1)
                ind = ind+1;
            end
        end
        ic = mod(ind-1,Ncolors)+1;
        for jj=1:sdip.n_dip
            im = mod(jj-1,Nmarker)+1;
            
            loc_pl = sdip.loc{tabl_seed_dip(ii,1)}(:,jj);
            js = sdip.j{tabl_seed_dip(ii,1)}(jj*3+l3,sdip.Mtb);
            vloc = sdip.cov_loc{tabl_seed_dip(ii,1)}(jj*3+l3,jj*3+l3);
            js_m = js_m+js;
            dip_h(:,ii+(jj-1)*pi_dip(1)) = ...
                add1dip(loc_pl,js/Mn_j*20,vloc, ...
                        marker{im},colors{ic},st.vols{1}.ax,Fig,st.bb);
        end
    end
end
st.vols{1}.sdip.ax = dip_h;

% Display dipoles orientation and strength
js_m = js_m/size(tabl_seed_dip,1);
[th,phi,Ijs_m] = cart2sph(js_m(1),js_m(2),js_m(3));
Njs_m = round(js_m'/Ijs_m*100)/100;
Angle = round([th phi]*1800/pi)/10;

set(sdip.hdl.hor1,'String',[num2str(Njs_m(1)),' ',num2str(Njs_m(2)), ...
        ' ',num2str(Njs_m(3))]);
set(sdip.hdl.hor2,'String',[num2str(Angle(1)),'� ',num2str(Angle(2)),'�']);
set(sdip.hdl.int,'String',Ijs_m);

% Change the colour of toggle button of dipoles actually displayed
for ii=tabl_seed_dip(:,1)
    set(sdip.hdl.hseed(ii),'BackgroundColor',[.7 1 .7]);
end


%==========================================================================
case 'clearall'
%==========================================================================
% Clears all dipoles, and reset the toggle buttons
%--------------------------------------------------------------------------

if isfield(st.vols{1},'sdip')
    sdip = st.vols{1}.sdip;
else
    error('I can''t find sdip structure');
end

disp('Clears all dipoles')
spm_eeg_inv_vbecd_disp('ClearDip');
for ii=1:st.vols{1}.sdip.n_seeds
    if sdip.exitflag(ii)==1
        set(st.vols{1}.sdip.hdl.hseed(ii),'Value',0);
    end
end
set(st.vols{1}.sdip.hdl.hdip,'Value',1);


%==========================================================================
case 'chgseed'
%==========================================================================
% Changes the seeds displayed
%--------------------------------------------------------------------------

% disp('Change seed')

sdip = st.vols{1}.sdip;
if isfield(sdip,'tabl_seed_dip')
    prev_seeds = p_seed(sdip.tabl_seed_dip);
else
    prev_seeds = [];
end

l_seed = zeros(sdip.n_seeds,1);
for ii=1:sdip.n_seeds
    if sdip.exitflag(ii)==1
        l_seed(ii) = get(sdip.hdl.hseed(ii),'Value');
    end
end
l_seed = find(l_seed);

% Modify the list of seeds displayed
if isempty(l_seed)
    % Nothing left displayed
    i_seed=[];
elseif isempty(prev_seeds)
    % Just one dipole added, nothing before
    i_seed=l_seed;
elseif length(prev_seeds)>length(l_seed)
    % One seed removed
    i_seed = prev_seeds;
    for ii=1:length(l_seed)
        p = find(prev_seeds==l_seed(ii));
        if ~isempty(p)
            prev_seeds(p) = [];
        end % prev_seeds is left with the index of the one removed
    end
    i_seed(i_seed==prev_seeds) = [];
    % Remove the dipole & change the button colour
    spm_eeg_inv_vbecd_disp('ClearDip',prev_seeds);
    set(sdip.hdl.hseed(prev_seeds),'BackgroundColor',[.7 .7 .7]);
else
    % One dipole added
    i_seed = prev_seeds;
    for ii=1:length(prev_seeds)
        p = find(prev_seeds(ii)==l_seed);
        if ~isempty(p)
            l_seed(p) = [];
        end % l_seed is left with the index of the one added
    end
    i_seed = [i_seed ; l_seed];
end

i_dip = get(sdip.hdl.hdip,'Value');
spm_eeg_inv_vbecd_disp('ClearDip');
if ~isempty(i_seed)
    spm_eeg_inv_vbecd_disp('DrawDip',i_seed,i_dip);
end


%==========================================================================
case 'chgdip'
%==========================================================================
% Changes the dipole index for the first seed displayed
%--------------------------------------------------------------------------

disp('Change dipole')
sdip = st.vols{1}.sdip;

i_dip = get(sdip.hdl.hdip,'Value');
if isfield(sdip,'tabl_seed_dip')
    i_seed = p_seed(sdip.tabl_seed_dip);
else
    i_seed = [];
end

if ~isempty(i_seed)
    spm_eeg_inv_vbecd_disp('ClearDip')
    spm_eeg_inv_vbecd_disp('DrawDip',i_seed,i_dip);
end


%==========================================================================
case 'cleardip'
%==========================================================================
% FORMAT spm_eeg_inv_vbecd_disp('ClearDip',seed_i)
% e.g. spm_eeg_inv_vbecd_disp('ClearDip')
%       clears all displayed dipoles
% e.g. spm_eeg_inv_vbecd_disp('ClearDip',1)
%       clears the first dipole displayed
%--------------------------------------------------------------------------

if nargin>2
    seed_i = varargin{1};
else
    seed_i = 0;
end
if isfield(st.vols{1},'sdip')
        sdip = st.vols{1}.sdip;
    else
        return; % I don't do anything, as I can't find sdip strucure
end

if isfield(sdip,'ax')
    Nax = size(sdip.ax,2);
    else
        return; % I don't do anything, as I can't find axes info
end

if seed_i==0 % removes everything
    for ii=1:Nax
        for jj=1:9
            delete(sdip.ax(jj,ii));
        end
    end
    for ii=sdip.tabl_seed_dip(:,1)
        set(sdip.hdl.hseed(ii),'BackgroundColor',[.7 .7 .7]);
    end
    sdip = rmfield(sdip,'tabl_seed_dip');
    sdip = rmfield(sdip,'ax');    
elseif seed_i<=Nax % remove one seed only
    l_seed = find(sdip.tabl_seed_dip(:,1)==seed_i);
    for ii=l_seed
        for jj=1:9
            delete(sdip.ax(jj,ii));
        end
    end
    sdip.ax(:,l_seed) = [];
    sdip.tabl_seed_dip(l_seed,:) = [];
else
    error('Trying to clear unspecified dipole');
end
st.vols{1}.sdip = sdip;


%==========================================================================
case 'redrawdip'
%==========================================================================
% spm_eeg_inv_vbecd_disp('RedrawDip')
% redraw everything, useful when zooming into image
%--------------------------------------------------------------------------

% spm_eeg_inv_vbecd_disp('ClearDip')
% spm_eeg_inv_vbecd_disp('ChgDip')

% disp('Change dipole')
sdip = st.vols{1}.sdip;

i_dip = get(sdip.hdl.hdip,'Value');
if isfield(sdip,'tabl_seed_dip')
    i_seed = p_seed(sdip.tabl_seed_dip);
else
    i_seed = [];
end

if ~isempty(i_seed)
    spm_eeg_inv_vbecd_disp('ClearDip')
    spm_eeg_inv_vbecd_disp('DrawDip',i_seed,i_dip);
end


%==========================================================================
otherwise
%==========================================================================
    warning('Unknown action string');
end
% warning(sw);
return


%==========================================================================
% dh = add1dip(loc,js,vloc,mark,col,ax,Fig,bb)
%==========================================================================
function dh = add1dip(loc,js,vloc,mark,col,ax,Fig,bb)
% Plots the dipoles on the 3 views, with an error ellipse for location
% Then returns the handle to the plots

global st
is = inv(st.Space);
loc = is(1:3,1:3)*loc(:) + is(1:3,4);
% taking into account the zooming/scaling only for the location
% NOT for the dipole's amplitude.
% Amplitude plotting is quite arbitrary anyway and up to some scaling
% defined for better viewing...

loc(1,:) = loc(1,:) - bb(1,1)+1;
loc(2,:) = loc(2,:) - bb(1,2)+1;
loc(3,:) = loc(3,:) - bb(1,3)+1;
% +1 added to be like John's orthview code

% prepare error ellipse
vloc = is(1:3,1:3)*vloc*is(1:3,1:3);
[V,E] = eig(vloc);
VE = V*diag(sqrt(diag(E))); % use std
% VE = V*E;   % or use variance ???

dh = zeros(9,1);
figure(Fig)
% Transverse slice, # 1
%----------------------
set(Fig,'CurrentAxes',ax{1}.ax)
set(ax{1}.ax,'NextPlot','add')
dh(1) = plot(loc(1),loc(2),[mark,col],'LineWidth',1);
dh(2) = plot(loc(1)+[0 js(1)],loc(2)+[0 js(2)],col,'LineWidth',2);
% add error ellipse
[uu,ss,vv] = svd(VE([1 2],:));
[phi] = cart2pol(uu(1,1),uu(2,1));
e = diag(ss);
t = (-1:.02:1)*pi;
x = e(1)*cos(t)*cos(phi)-e(2)*sin(t)*sin(phi)+loc(1);
y = e(2)*sin(t)*cos(phi)+e(1)*cos(t)*sin(phi)+loc(2);
dh(3) = plot(x,y,[':',col],'LineWidth',.5);
set(ax{1}.ax,'NextPlot','replace')

% Coronal slice, # 2
%----------------------
set(Fig,'CurrentAxes',ax{2}.ax)
set(ax{2}.ax,'NextPlot','add')
dh(4) = plot(loc(1),loc(3),[mark,col],'LineWidth',1);
dh(5) = plot(loc(1)+[0 js(1)],loc(3)+[0 js(3)],col,'LineWidth',2);
% add error ellipse
[uu,ss,vv] = svd(VE([1 3],:));
[phi] = cart2pol(uu(1,1),uu(2,1));
e = diag(ss);
t = (-1:.02:1)*pi;
x = e(1)*cos(t)*cos(phi)-e(2)*sin(t)*sin(phi)+loc(1);
y = e(2)*sin(t)*cos(phi)+e(1)*cos(t)*sin(phi)+loc(3);
dh(6) = plot(x,y,[':',col],'LineWidth',.5);
set(ax{2}.ax,'NextPlot','replace')

% Sagital slice, # 3
%----------------------
set(Fig,'CurrentAxes',ax{3}.ax)
set(ax{3}.ax,'NextPlot','add')
% dh(5) = plot(dim(2)-loc(2),loc(3),[mark,col],'LineWidth',2);
% dh(6) = plot(dim(2)-loc(2)+[0 -js(2)],loc(3)+[0 js(3)],col,'LineWidth',2);
dh(7) = plot(bb(2,2)-bb(1,2)-loc(2),loc(3),[mark,col],'LineWidth',1);
dh(8) = plot(bb(2,2)-bb(1,2)-loc(2)+[0 -js(2)],loc(3)+[0 js(3)],col,'LineWidth',2);
% add error ellipse
[uu,ss,vv] = svd(VE([2 3],:));
[phi] = cart2pol(uu(1,1),uu(2,1));
e = diag(ss);
t = (-1:.02:1)*pi;
x = -(e(1)*cos(t)*cos(phi)-e(2)*sin(t)*sin(phi))+bb(2,2)-bb(1,2)-loc(2);
y = e(2)*sin(t)*cos(phi)+e(1)*cos(t)*sin(phi)+loc(3);
dh(9) = plot(x,y,[':',col],'LineWidth',.5);
set(ax{3}.ax,'NextPlot','replace')

return

%==========================================================================
% pr_seed = p_seed(tabl_seed_dip)
%==========================================================================
function pr_seed = p_seed(tabl_seed_dip)
% Gets the list of seeds used in the previous display

ls = sort(tabl_seed_dip(:,1));
if length(ls)==1
    pr_seed = ls;
else
    pr_seed = ls([find(diff(ls)) ; length(ls)]);
end

%


% OLD STUFF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use it with arguments or not:
% - spm_eeg_inv_vbecd_disp('Init')
%        The routine asks for the dipoles file and image to display
% - spm_eeg_inv_vbecd_disp('Init',sdip)
%        The routine will use the avg152T1 canonical image 
% - spm_eeg_inv_vbecd_disp('Init',sdip,P)
%        The routines dispays the dipoles on image P.
%
% If multiple seeds have been used, you can select the seeds to display 
% by pressing their index. 
% Given that the sources could have different locations, the slices
% displayed will be the 3D view at the *average* or *mean* locations of
% selected sources.
% If more than 1 dipole was fitted at a time, then selection of source 1 
% to N is possible through the pull-down selector.
%
% The location of the source/cut is displayed in mm and voxel, as well as
% the underlying image intensity at that location.
% The cross hair position can be hidden by clicking on its button.
%
% Nota_1: If the cross hair is manually moved by clicking in the image or
%       changing its coordinates, the dipole displayed will NOT be at
%       the right displayed location. That's something that needs to be improved...
%
% Nota_2: Some seeds may have not converged within the limits fixed,
%       these dipoles are not displayed...
%
% Fields needed in sdip structure to plot on an image:
%       + n_seeds: nr of seeds set used, i.e. nr of solutions calculated
%       + n_dip: nr of fitted dipoles on the EEG time series
%       + loc: location of fitted dipoles, cell{1,n_seeds}(3 x n_dip)
%               remember that loc is fixed over the time window.
%       + j: sources amplitude over the time window, 
%            cell{1,n_seeds}(3*n_dip x Ntimebins)
%       + Mtb: index of maximum power in EEG time series used


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % First point to consider
% loc_mm = sdip.loc{i_seed(1)}(:,i_dip);
% 
% % PLace the underlying image at right cuts
% spm_orthviews('Reposition',loc_mm);
% % spm_orthviews('Reposition',loc_vx);
% % spm_orthviews('Xhairs','off')
% 
% % if i_seed = set, Are there other dipoles close enough ?
% tabl_seed_dip=[i_seed(1) i_dip]; % table summarising which set & dip to use.
% if length(i_seed)>1
%   unit = ones(1,sdip.n_dip);
%   for ii = i_seed(2:end)'
%         d2 = sqrt(sum((sdip.loc{ii}-loc_mm*unit).^2));
%         l_cl = find(d2<=lim_cl);
%         if ~isempty(l_cl)
%             for jj=l_cl
%                 tabl_seed_dip = [tabl_seed_dip ; [ii jj]];
%             end
%         end
%   end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get(sdip.hdl.hseed(1),'Value')
% for ii=1:sdip.n_seeds, delete(hseed(ii)); end
% h1 = uicontrol(Fig,'Style','togglebutton','Position',[600 25 10 10].*WS)
% h2 = uicontrol(Fig,'Style','togglebutton','Position',[620 100 20 20].*WS,'String','1')
% h2 = uicontrol(Fig,'Style','checkbox','Position',[600 100 10 10].*WS)
% h3 = uicontrol(Fig,'Style','radiobutton','Position',[600 150 20 20].*WS)
% h4 = uicontrol(Fig,'Style','radiobutton','Position',[700 150 20 20].*WS)
% delete(h2),delete(h3),delete(h4),
% delete(hdip)
