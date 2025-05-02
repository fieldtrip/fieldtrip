%  HBF_INTERPOLATETFULLTOELECTRODES extracts the transfer matrix for
%  potential at given electrode locations (on the outer boundary only)
% 
%  function Tphi_elecs=HBF_INTERPOLATETFULLTOELECTRODES(Tfull,bmeshes,elecs);
%    Tfull:    full bf BEM transfer matrix
%    bmeshes:  hbf bmeshes struct
%    elecs:    hbf electrodes struct
% 
%    Tphi_elecs: transfer matrix for computing potential at electrode
%    locations
% 
%    v200924 Matti Stenroos
%