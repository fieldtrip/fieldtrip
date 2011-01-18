function spm_normalise_disp(matname,VF)
% Display results of spatial normalisation
% FORMAT spm_normalise_disp(matname)
% matname - name of sn3d.mat file
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_normalise_disp.m 1143 2008-02-07 19:33:33Z spm $


fg = spm_figure('FindWin','Graphics');
if isempty(fg), return; end;

if nargin<1, matname = spm_select(1,'.*_sn.mat$','Select parameter file'); end;

if ischar(matname),
    t = load(deblank(matname));
else %assume it is a structure
    t = matname;
end;

if nargin<2, VF = t.VF(1); end;

Q = t.VG(1).mat*inv(t.Affine)/VF.mat;

spm_figure('Clear','Graphics');
ax=axes('Position',[0.1 0.51 0.8 0.45],'Visible','off','Parent',fg);
text(0,0.90, 'Spatial Normalisation','FontSize',16,'FontWeight','Bold',...
    'Interpreter','none','Parent',ax);
text(0,0.75, [ 'Image     : ' VF.fname],'FontSize',12,'FontWeight','Bold',...
    'Interpreter','none','Parent',ax);
%text(0,0.7, [ 'Parameters : ' spm_str_manip(matname,'sd')],'FontSize',12,...
%   'Interpreter','none','Parent',ax);

%str = 'no flipping';
%if det(t.Affine(1:3,1:3))<0, str = 'image flipped'; end;
text(0,0.6, 'Linear {affine} component','FontWeight','Bold',...
    'Interpreter','none','Parent',ax);
text(0,0.55, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),...
    'Interpreter','none','Parent',ax);
text(0,0.50, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),...
    'Interpreter','none','Parent',ax);
text(0,0.45, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),...
    'Interpreter','none','Parent',ax);

d = [size(t.Tr) 1 1 1];
d = d(1:3);

if prod(d)>1 && isfinite(t.flags.reg),
    text(0,0.35, sprintf('%d nonlinear iterations',t.flags.nits),...
        'Interpreter','none','Parent',ax);
    text(0,0.30, sprintf('%d x %d x %d basis functions',d),...
        'Interpreter','none','Parent',ax);
else
    text(0,0.35, 'No nonlinear components',...
        'Interpreter','none','Parent',ax);
end;

spm_orthviews('Reset');
spm_orthviews('Image',t.VG(1).fname,[0.01 0.1 .48 .6]);
VN = spm_write_sn(VF,matname);
h2 = spm_orthviews('Image',VN,[.51 0.1 .48 .6]);
spm_orthviews('Space',h2);
spm_print;
drawnow;
return;

