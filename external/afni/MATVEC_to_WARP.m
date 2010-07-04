Wd = function (matvec)
%take a 12 param matvec, typically WARPDRIVE_MATVEC_INV_000000 and
%change it to a WARP_DATA
%Used to test transformations needed to make AFNI use 3dWarpDrive's output
%as a TLRC transform

if (nargin == 0) matvec = [];

if (isempty(matvec)),
   WARNING: Using preset matvecs!\n\n');
   %test case, matvecs from DemoSubj_SurfVol_Alnd_Exp_at2+tlrc
   if (0),
      matvec = [  0.9709206    -0.04199445     0.01165925     -0.7014965     
               0.07297896    0.9860933    -0.08915191     -0.2648731    
               0.005342953    -0.02306779  0.9538299      -6.124281
               ]; 
   else
      %that would be WARPDRIVE_MATVEC_INV_000000
      %the proper one to use for WARP_DATA, the one that takes you 
      %from +orig to +tlrc
      matvec = [  1.026725     0.04352641   -0.008482005      0.6798272    
                  -0.07667368  1.013075     0.09562658      0.8001939   
                  -0.007605589   0.02425676 1.050765       6.436271
               ]; 
   end
else
   if (length(matvec) ~= 12),
      fprintf(2,'Error: matvec must be 12 elements long\n');
   end
end

Mfor = [matvec; 0 0 0 1]; %the forward transform
Mbac = inv(Mfor); %the inverse transform

%begin transformation
mfor = Mfor(1:3, 1:3);
mbac = Mbac(1:3, 1:3);
bvec = -Mfor(1:3,4); %Don't ask about the minus sign, see README.attributes:
                     %Under WARP_DATA, The forward transformation is  [x_map] = [mfor] [x_in]  - [bvec];
                     %Under Registration Attributes: [xyz_out] = [mat] ([xyz_in]-[xyz_cen]) + [vec] + [xyz_cen]
svec = -Mbac(1:3,4);
bot = [-80 ; -80 ; -65]; %good enough for TT box
top = [80 ; 110 ; 85]; 

Wd = [reshape(mfor',1,9) reshape(mbac',1,9) reshape(bvec,1,3) reshape(svec,1,3) [bot'] [top']];

fprintf(1,'Wd = \n\n'); 
fprintf(1,'       %g   %g   %g   %g   %g   %g\n',Wd)
fprintf(1,'\n\n'); 
