function [transfer] = sb_transfer(vol,sens)

% SB_TRANSFER
% 
% INPUT:
% vol.wf.nd,vol.wf.el: nodes and elements of head grid (nodes in mm)
% vol.wf.field: vector assigning a conductivity to teach element (in S/mm)
% elc: vector containing the electrode positions (in mm)
% 
% OUTPUT:
% vol.transfer: transfer matrix
%
% $Id$

%--------------------------------------------------------------------------
%TODO: tissuecond: vector with the conductivities of the respective labels,
%right now not used
%--------------------------------------------------------------------------
%TODO: automated setting of conductivities?
%disp('Setting conductivities...')
%vol.wf.field = sb_set_cond(vol.wf.labels,tissuecond); %gucken, wo tissuecond herkommt...
%--------------------------------------------------------------------------
%TODO: this probably needs to be expanded so that the node is a surface
%node if no ECoG-flag is set - solution: find a surface first...
%--------------------------------------------------------------------------

disp('Find electrode positions...')
vol.elecnodes = sb_find_elec(vol,sens);
%calculate transfermatrix
disp('Calculate transfer matrix...')
transfer = zeros(length(vol.elecnodes),size(vol.pos,1));
%--------------------------------------------------------------------------
%TODO: add possibility for parallel computation
%matlabpool local 2;
for i=2:length(vol.elecnodes)
    % NOTE: the loop starts at 2, which seems intentional, because the reference is the first electrode
    str = ['Electrode ',num2str(i),' of ',num2str(size(vol.elecnodes,1))];
    disp(str)
    vecb = zeros(size(vol.stiff,1),1);
    vecb(vol.elecnodes(i)) = 1;
    transfer(i,:) = sb_calc_vecx(vol.stiff,vecb,vol.elecnodes(1));
    clear vecb;
end
%matlabpool close;
end
