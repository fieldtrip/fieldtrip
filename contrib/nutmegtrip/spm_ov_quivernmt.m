function ret = spm_ov_quivernmt(varargin)
% This routine is a plugin to spm_orthviews for SPM12 / SPM8.
% Used by NutmegTrip / nmt_sourceoriplot.m to plot vector fields
% (e.g., orientations) on MRI
%
% For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% adapted from spm_ov_quiver.m disributed with SPM Tools / tbxDiffusion:
% https://sourceforge.net/projects/spmtools/



global st;

if isempty(st)
    error('quivernmt: This routine can only be called as a plugin for spm_orthviews!');
end;

if nargin < 2
    error('quivernmt: Wrong number of arguments. Usage: spm_orthviews(''quivernmt'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};
switch cmd
    case 'init'     %function addquiver(handle, Vqfnames, Vmaskfname, varargin)
        if nargin < 4
            error('spm_orthviews(''quivernmt'', ''init'',...): Not enough arguments');
        end;
        Vq = spm_vol(varargin{3});
        Vmask = spm_vol(varargin{4});
        if length(Vq) == 3
            st.vols{volhandle}.quivernmt = struct('qx',Vq(1),'qy',Vq(2),'qz',Vq(3), ...
                'mask',Vmask, 'fa',[], 'thresh', [.1 Inf], 'ls','b.', ...
                'qst',3,'ql',1,'qht',[],'qhc',[],'qhs',[], 'qw', .5);
        else
            error('spm_orthviews(''quivernmt'', ''init'',...): Please specify 3 images!');
        end;
        if nargin > 4
            if ~isempty(varargin{5})
                st.vols{volhandle}.quivernmt.fa = spm_vol(varargin{5});
            end;
        end;
        if nargin > 5
            if ~isempty(varargin{6})
                st.vols{volhandle}.quivernmt.thresh(1) = varargin{6};
            end;
        end;
        if nargin > 6
            if ~isempty(varargin{7})
                st.vols{volhandle}.quivernmt.thresh(2) = varargin{7};
            end;
        end;
        if nargin > 7
            if ischar(varargin{8})
                st.vols{volhandle}.quivernmt.ls = varargin{8};
            end;
        end;
        if nargin > 8
            if ~isempty(varargin{9})
                st.vols{volhandle}.quivernmt.qst = varargin{9};
            end;
        end;
        if isempty(st.vols{volhandle}.quivernmt.fa)
            st.vols{volhandle}.quivernmt.ql = 1;
        end;
        if nargin > 9
            if ~isempty(varargin{10})
                st.vols{volhandle}.quivernmt.ql = varargin{10};
            end;
        end;
        if nargin > 10
            if ~isempty(varargin{11})
                st.vols{volhandle}.quivernmt.qw = varargin{11};
            end;
        end;

    case 'redraw'
        TM0 = varargin{3};
        TD  = varargin{4};
        CM0 = varargin{5};
        CD  = varargin{6};
        SM0 = varargin{7};
        SD  = varargin{8};
        if isfield(st.vols{volhandle},'quivernmt')
            % need to delete old quiver lines before redrawing
            delete(st.vols{volhandle}.quivernmt.qht);
            delete(st.vols{volhandle}.quivernmt.qhc);
            delete(st.vols{volhandle}.quivernmt.qhs);

            qx  = st.vols{volhandle}.quivernmt.qx;
            qy  = st.vols{volhandle}.quivernmt.qy;
            qz  = st.vols{volhandle}.quivernmt.qz;
            mask = st.vols{volhandle}.quivernmt.mask;
            fa = st.vols{volhandle}.quivernmt.fa;
            thresh(1) = st.vols{volhandle}.quivernmt.thresh(1);
            thresh(2) = st.vols{volhandle}.quivernmt.thresh(2);
            ls = st.vols{volhandle}.quivernmt.ls;

            % step size for selection of locations
            prm = spm_imatrix(st.Space);
            qst = ceil(st.vols{volhandle}.quivernmt.qst/prm(7));
            qst = st.vols{1}.quivernmt.qx.mat(1,1); % get voxel/grid spacing from the transformation matrix
            ql = st.vols{volhandle}.quivernmt.ql; % scaling of arrow length
            qst1 = ceil(qst/2);

            Mx   = st.vols{volhandle}.premul*qx.mat;
            My   = st.vols{volhandle}.premul*qy.mat;
            Mz   = st.vols{volhandle}.premul*qz.mat;
            Mm   = st.vols{volhandle}.premul*mask.mat;
            if ~isempty(fa)
                fat = spm_slice_vol(fa,inv(TM0*(st.Space\Mx)),TD,0)';
            else
                fat = 1;
            end;
            rqt = cat(3, spm_slice_vol(qx,inv(TM0*(st.Space\Mx)),TD,0)', ...
                spm_slice_vol(qy,inv(TM0*(st.Space\My)),TD,0)', ...
                spm_slice_vol(qz,inv(TM0*(st.Space\Mz)),TD,0)');
            rqt = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(rqt,TD(1)*TD(2),3)';
            qxt = fat.*reshape(rqt(1,:)',TD(2),TD(1));
            qyt = fat.*reshape(rqt(2,:)',TD(2),TD(1));
            qzt = fat.*reshape(rqt(3,:)',TD(2),TD(1));

            maskt = spm_slice_vol(mask,inv(TM0*(st.Space\Mm)),TD,0)';
            xt = (1:TD(1))-.5;
            yt = (1:TD(2))-.5;
            zt = zeros(size(qxt));
            zt((maskt < thresh(1))|(maskt > thresh(2))) = NaN;
            zt((qxt == 0) & (qyt == 0) & (qzt == 0)) = NaN;

            if ~isempty(fa)
                fac = spm_slice_vol(fa,inv(CM0*(st.Space\Mx)),CD,0)';
            else
                fac = 1;
            end;
            rqc = cat(3, spm_slice_vol(qx,inv(CM0*(st.Space\Mx)),CD,0)', ...
                spm_slice_vol(qy,inv(CM0*(st.Space\My)),CD,0)', ...
                spm_slice_vol(qz,inv(CM0*(st.Space\Mz)),CD,0)');
            rqc = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(rqc,CD(1)*CD(2),3)';
            qxc = fac.*reshape(rqc(1,:)',CD(2),CD(1));
            qyc = fac.*reshape(rqc(2,:)',CD(2),CD(1));
            qzc = fac.*reshape(rqc(3,:)',CD(2),CD(1));

            maskc = spm_slice_vol(mask,inv(CM0*(st.Space\Mm)),CD,0)';
            xc = (1:CD(1))-.5;
            yc = (1:CD(2))-.5;
            zc = zeros(size(qxc));
            zc((maskc < thresh(1))|(maskc > thresh(2))) = NaN;
            zc((qxc == 0) & (qyc == 0) & (qzc == 0)) = NaN;

            if ~isempty(fa)
                fas = spm_slice_vol(fa,inv(SM0*(st.Space\Mx)),SD,0)';
            else
                fas = 1;
            end;
            rqs = cat(3, spm_slice_vol(qx,inv(SM0*(st.Space\Mx)),SD,0)', ...
                spm_slice_vol(qy,inv(SM0*(st.Space\My)),SD,0)', ...
                spm_slice_vol(qz,inv(SM0*(st.Space\Mz)),SD,0)');
            rqs = st.Space(1:3,1:3)*st.vols{volhandle}.premul(1:3,1:3)*reshape(rqs,SD(1)*SD(2),3)';
            qxs = fas.*reshape(rqs(1,:)',SD(2),SD(1));
            qys = fas.*reshape(rqs(2,:)',SD(2),SD(1));
            qzs = fas.*reshape(rqs(3,:)',SD(2),SD(1));

            masks = spm_slice_vol(mask,inv(SM0*(st.Space\Mm)),SD,0)';
            xs = (1:SD(1))-.5;
            ys = (1:SD(2))-.5;
            zs = zeros(size(qxs));
            zs((masks < thresh(1))|(masks > thresh(2))) = NaN;
            zs((qxs == 0) & (qys == 0) & (qzs == 0)) = NaN;

            % check for availability of "centered" quiver function
            if exist('nmt_tbxdti_quiver3.m','file')
                quiverfun = @nmt_tbxdti_quiver3;
            else % fallback
                quiverfun = @quiver3;
                warning('Function "tbxdti_quiver3.m" not found!\n Using standard Matlab routine "quiver3.m" instead.\n This may not give  nice quiver plots.');
            end;
            if spm_flip_analyze_images
                flipx = -1;
            else
                flipx = 1;
            end;

            % transversal - plot (x y z)
            np = get(st.vols{volhandle}.ax{1}.ax,'NextPlot');
            set(st.vols{volhandle}.ax{1}.ax,'NextPlot','add');
            axes(st.vols{volhandle}.ax{1}.ax);
            x=xt(qst1:qst:end);
            y=yt(qst1:qst:end);
            z=zt(qst1:qst:end,qst1:qst:end);
            u=flipx*qxt(qst1:qst:end,qst1:qst:end);
            v=      qyt(qst1:qst:end,qst1:qst:end);
            w=      qzt(qst1:qst:end,qst1:qst:end);
            xtt=repmat(x,size(z,1),1);
            ytt=repmat(y',1,size(z,2));
            st.vols{volhandle}.quivernmt.qht = feval(quiverfun,xtt,ytt,z,u,v,w,ql,ls);
            set(st.vols{volhandle}.ax{1}.ax,'NextPlot',np);
            set(st.vols{volhandle}.quivernmt.qht, ...
                'Parent',st.vols{volhandle}.ax{1}.ax, 'HitTest','off', ...
                'Linewidth',st.vols{volhandle}.quivernmt.qw );

            % coronal - plot (x z y)
            np = get(st.vols{volhandle}.ax{2}.ax,'NextPlot');
            set(st.vols{volhandle}.ax{2}.ax,'NextPlot','add');
            axes(st.vols{volhandle}.ax{2}.ax);
            x=xc(qst1:qst:end);
            y=yc(qst1:qst:end);
            z=zc(qst1:qst:end,qst1:qst:end);
            u=flipx*qxc(qst1:qst:end,qst1:qst:end);
            v=      qzc(qst1:qst:end,qst1:qst:end);
            w=      qyc(qst1:qst:end,qst1:qst:end);
            xtt=repmat(x,size(z,1),1);
            ytt=repmat(y',1,size(z,2));
            st.vols{volhandle}.quivernmt.qhc = feval(quiverfun,xtt,ytt,z,u,v,w,ql,ls);
            set(st.vols{volhandle}.ax{2}.ax,'NextPlot',np);
            set(st.vols{volhandle}.quivernmt.qhc, ...
                'Parent',st.vols{volhandle}.ax{2}.ax, 'HitTest','off', ...
                'Linewidth',st.vols{volhandle}.quivernmt.qw );

            % sagittal - plot (-y z x)
            np = get(st.vols{volhandle}.ax{3}.ax,'NextPlot');
            set(st.vols{volhandle}.ax{3}.ax,'NextPlot','add');
            axes(st.vols{volhandle}.ax{3}.ax);
            x=xs(qst1:qst:end);
            y=ys(qst1:qst:end);
            z=zs(qst1:qst:end,qst1:qst:end);
            u=     -qys(qst1:qst:end,qst1:qst:end);
            v=      qzs(qst1:qst:end,qst1:qst:end);
            w=flipx*qxs(qst1:qst:end,qst1:qst:end);
            xtt=repmat(x,size(z,1),1);
            ytt=repmat(y',1,size(z,2));
            st.vols{volhandle}.quivernmt.qhs = feval(quiverfun,xtt,ytt,z,u,v,w,ql,ls);
            set(st.vols{volhandle}.ax{3}.ax,'NextPlot',np);
            set(st.vols{volhandle}.quivernmt.qhs, ...
                'Parent',st.vols{volhandle}.ax{3}.ax, 'HitTest','off', ...
                'Linewidth',st.vols{volhandle}.quivernmt.qw );
        end; %quiver

    case 'delete'
        if isfield(st.vols{volhandle},'quivernmt'),
            delete(st.vols{volhandle}.quivernmt.qht);
            delete(st.vols{volhandle}.quivernmt.qhc);
            delete(st.vols{volhandle}.quivernmt.qhs);
            st.vols{volhandle} = rmfield(st.vols{volhandle},'quivernmt');
        end;
        %-------------------------------------------------------------------------
        % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, 'Label', 'Quiver');
        item1 = uimenu(item0, 'Label', 'Add', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_init'', ', ...
            num2str(volhandle), ');'], 'Tag', ['QUIVER_0_', num2str(volhandle)]);
        item2 = uimenu(item0, 'Label', 'Properties', ...
            'Visible', 'off', 'Tag', ['QUIVER_1_', num2str(volhandle)]);
        item2_1 = uimenu(item2, 'Label', 'Mask threshold', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_edit'',', ...
            num2str(volhandle), ',''thresh'');']);
        item2_2 = uimenu(item2, 'Label', 'Linestyle', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_edit'',', ...
            num2str(volhandle), ',''ls'');']);
        item2_3 = uimenu(item2, 'Label', 'Quiver distance', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_edit'',', ...
            num2str(volhandle), ',''qst'');']);
        item2_4 = uimenu(item2, 'Label', 'Quiver length', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_edit'',', ...
            num2str(volhandle), ',''ql'');']);
        item2_5 = uimenu(item2, 'Label', 'Linewidth', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_edit'',', ...
            num2str(volhandle), ',''qw'');']);
        item3 = uimenu(item0, 'Label', 'Remove', 'Callback', ...
            ['feval(''spm_ov_quivernmt'',''context_delete'', ', ...
            num2str(volhandle), ');'], 'Visible', 'off', ...
            'Tag', ['QUIVER_1_', num2str(volhandle)]);

    case 'context_init'
        Finter = spm_figure('FindWin', 'Interactive');
        spm_input('!DeleteInputObj',Finter);
        [Vqfnames, sts] = spm_select(3, 'image',...
            'Components of 1st eigenvector',[], ...
            pwd, 'evec1.*');
        if ~sts
            return;
        end;
        [Vmaskfname, sts] = spm_select(1,'image', 'Mask image');
        if ~sts
            return;
        end;
        Vfafname = spm_select(Inf,'image', 'Fractional anisotropy image', [], ...
            pwd, 'fa.*');
        feval('spm_ov_quivernmt','init',volhandle,Vqfnames,Vmaskfname,Vfafname);
        obj = findobj(0, 'Tag',  ['QUIVER_1_', num2str(volhandle)]);
        set(obj, 'Visible', 'on');
        obj = findobj(0, 'Tag',  ['QUIVER_0_', num2str(volhandle)]);
        set(obj, 'Visible', 'off');
        spm_orthviews('redraw');

    case 'context_edit'
        Finter = spm_figure('FindWin', 'Interactive');
        spm_input('!DeleteInputObj',Finter);
        switch varargin{3}
            case 'thresh'
                in = spm_input('Mask threshold {min max}','!+1','e', ...
                    num2str(st.vols{volhandle}.quivernmt.thresh), [1 2]);
            case 'ls'
                in = spm_input('Line style','!+1','s', ...
                    st.vols{volhandle}.quivernmt.ls);
            case 'qst'
                in = spm_input('Quiver distance','!+1','e', ...
                    num2str(st.vols{volhandle}.quivernmt.qst), 1);
            case 'ql'
                in = spm_input('Quiver length','!+1','e', ...
                    num2str(st.vols{volhandle}.quivernmt.ql), 1);
            case 'qw'
                in = spm_input('Linewidth','!+1','e', ...
                    num2str(st.vols{volhandle}.quivernmt.qw), 1);
        end;
        spm_input('!DeleteInputObj',Finter);
        st.vols{volhandle}.quivernmt.(varargin{3}) = in;
        spm_orthviews('redraw');

    case 'context_delete'
        feval('spm_ov_quivernmt','delete',volhandle);
        obj = findobj(0, 'Tag',  ['QUIVER_1_', num2str(volhandle)]);
        set(obj, 'Visible', 'off');
        obj = findobj(0, 'Tag',  ['QUIVER_0_', num2str(volhandle)]);
        set(obj, 'Visible', 'on');

    otherwise

        fprintf('spm_orthviews(''quiver'',...): Unknown action string %s', cmd);
end;
