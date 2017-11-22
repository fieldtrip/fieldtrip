function failed_fairsurface

% MEM 1500mb
% WALLTIME 00:10:00

ft_defaults

hdmfile  = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.shape');
shape = ft_read_headshape(hdmfile);

tmpcfg             = [];
tmpcfg.numvertices = [];
tmpcfg.headshape   = shape;
geometry1=ft_prepare_mesh(tmpcfg);
% output is
% geometry1 = 
%     pnt: [21285x3 double]
%     tri: [42566x3 double]

tmpcfg             = [];
tmpcfg.numvertices = 3000; % default set in ft_prepare_headmodel
tmpcfg.headshape   = shape;
geometry2=ft_prepare_mesh(tmpcfg);
% output is 
% geometry2 = 
%     pnt: [2634x3 double]
%     tri: [5204x3 double]
% but note that several geometry2.pnt are NaN    
if length(find(isnan(geometry2.pnt(:))))>0
  error
end
 
