function [transfer vol] = sb_transfer(vol,elc,tissuecond);
%INPUT:
%vol.wf.nd,vol.wf.el: nodes and elements of head grid (nodes in mm)
%vol.wf.field: vector assigning a conductivity to teach element (in S/mm)
%elc: vector containing the electrode positions (in mm)
%OUTPUT: vol.transfer: transfer matrix
%--------------------------------------------------------------------------
%TODO: tissuecond: vector with the conductivities of the respective labels,
%right now not used
%TODO: automated setting of conductivities?
%disp('Setting conductivities...')
%vol.wf.field = sb_set_cond(vol.wf.labels,tissuecond); %gucken, wo tissuecond herkommt...
%--------------------------------------------------------------------------
%stiffnes matrix
disp('Calculate stiffnes matrix...')
stiff = sb_calc_stiff(vol.wf.nd,vol.wf.el,vol.wf.field);
%find nodes corresponding to electrode positions
%--------------------------------------------------------------------------
%TODO: this probably needs to be expanded so that the node is a surface
%node if no ECoG-flag is set - solution: find a surface first...
disp('Find electrode positions...')
vol.wf.diri = sb_find_elec(vol,elc);
%calculate transfermatrix
disp('Calculate transfer matrix...')
transfer = zeros(size(elc,1),size(vol.wf.nd,1));
%--------------------------------------------------------------------------
%TODO: add possibility for parallel computation
%matlabpool local 2;
for i=2:length(elc)
    str = ['Electrode ',num2str(i),' of ',num2str(size(elc,1))];
    disp(str)
    vecb = zeros(size(stiff,1),1);
    vecb(vol.wf.diri(i)) = 1;
    transfer(i,:) = sb_calc_vecx(stiff,vecb,vol.wf.diri(1));
    clear vecb;
end
%matlabpool close;
%vol.transfer = transfer;
%clear transfer;
end
