% DeOblique
% a script that figures out the rotation to apply to
% cardinal-plane-acquired data to bring it into alignment with oblique functional data
%
% Usage:
% The script will prompt you for an Ifile filename, say: may0603mj/I.610
% and will output the 3drotate command needed to apply to the
% High resolution anatomical volume to bring it into alignment with
% the Ifile (may0603mj/I.610)
%
% Hit enter to stop the loop. The results are logged in a file called
% DeOblique.txt
%
%This script, with its silly name, was never meant to be distributed,
%it was hastily written and meant for my own consumption only.
%It has been only tested for slices slightly off the axial plane.
%
%If you are using this script for some reason, use it with
%caution. Verify your results with AFNI at the end.
%
%saadz@mail.nih.gov

FuncName = 'DeOblique';

%load image header
%RAS here means LPI

%fname = sprintf('/home/ziad/tmp/pm040803/5to10/023/I.138');
%fname = sprintf('/home/ziad/tmp/pm040803/6to10/024/I.138');
%fname = sprintf('/home/ziad/tmp/pm040803/spgr/005/I.073');



fname = input ('Enter I filename (or test or dircos, or twovec): ', 's');
while (~isempty(fname)),
   if (strcmp(fname,'test') == 1),
      tlhc1 = [83 -123 32];
      trhc1 = [-83 -123 32];
      brhc1 = [-83 -55 22];
      ctr1  = [ 0   0 27];
      Tri1.XYZ = [tlhc1; trhc1; brhc1];
      [err,Eq1] = Plane_Equation (Tri1);
   elseif (strcmp(fname,'dircos') == 1),
      d1 = []; while (length(d1) ~= 3), d1 = input ('Enter 1st direction cosines [1 0 0]: '); end
      d2 = []; while (length(d2) ~= 3), d2 = input ('Enter 2nd direction cosines [0 0.88 0.44]: '); end
      d1 = acos(d1); d2 = acos(d2);
      dpl = cross(d1, d2); % plane normal
      ctr1 = [0 0 0];
      Eq1 = [dpl 0];
      Tri1 = [];
   elseif (strcmp(fname,'twovec') == 1),
      d1 = []; while (length(d1) ~= 3), d1 = input ('Enter 1st vector [1 0 0]: '); end
      d2 = []; while (length(d2) ~= 3), d2 = input ('Enter 2nd vector [0 1 0]: '); end
      dpl = cross(d1, d2); % plane normal
      ctr1 = [0 0 0];
      Eq1 = [dpl 0];
      Tri1 = [];
   elseif (strcmp(fname, 'scan_6') == 1),
      tlhc1 = [125.994 	145.952 	-25.001];
      trhc1 = [-113.867 152.413  -29.908];
      brhc1 = [-120.694 -86.730  -10.843];
      ctr1  = [2.650 	29.611 	-17.922];
      Tri1.XYZ = [tlhc1; trhc1; brhc1];
      [err,Eq1] = Plane_Equation (Tri1);
    else
      [su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(fname);

      %get three points forming corners of plane and the center of the plane
      tlhc1 = [im_hdr.tlhc_R im_hdr.tlhc_A im_hdr.tlhc_S];
      trhc1 = [im_hdr.trhc_R im_hdr.trhc_A im_hdr.trhc_S];
      brhc1 = [im_hdr.brhc_R im_hdr.brhc_A im_hdr.brhc_S];
      ctr1 = [im_hdr.ctr_R im_hdr.ctr_A im_hdr.ctr_S];
      Tri1.XYZ = [tlhc1; trhc1; brhc1];
      [err,Eq1] = Plane_Equation (Tri1);
   end

   %an idea of the displacement at the edges, typically the minimum of d
   if (~isempty(Tri1)),
      d = 2 .* (tlhc1 - ctr1);
   else
      d = -1;
   end

   %form the equation of the desired plane
   ctr2 = ctr1;

   for (car=1:1:3),
      if (car == 1) CardPlane = sprintf('Axial');
      elseif (car == 3) CardPlane = sprintf('Sagittal');
      else CardPlane = sprintf('Coronal');
      end
   %desired plane is one that is parallel to Axial plane and passing through ctr2
   switch (CardPlane(1)),
      case 'A',
         % Axial plane eqation is z - cnst = 0;
         Eq2 = [0 0 1 -ctr2(3)];
      case 'C',
         % Coronal plane equation is y - cnst = 0;
         Eq2 = [0 1 0 -ctr2(3)];
      case 'S',
         % Coronal plane equation is x - cnst = 0;
         Eq2 = [1 0 0 -ctr2(3)];
      otherwise,
         fprintf (1,'Error Bad Plane.\n');
   end

   % find the intersection line of the two planes
   u  = cross(Eq1(1:3),Eq2(1:3)); %the direction vector of the intersection line
   %should check that cross product is not ~ 0 (i.e. // planes)
   u = u ./ norm(u, 2);

   %This line must pass by ctr because I said so
   if (~isempty(Tri1)),
      somescale = abs(max(trhc1 - brhc1))./3;
   else
      somescale = 10;
   end
   Line = [ctr1; ctr1+somescale.*u];


   %Good, now get the angle
   CosAngle = dot(Eq1(1:3),Eq2(1:3)) ./ ( norm(Eq1(1:3), 2) .* norm(Eq2(1:3), 2));
   AngleRad = acos(CosAngle)
   AngleDeg = AngleRad .* 180 ./ pi

   %Show the results
   Opt.Fig = figure(car); clf

   subplot 211; cla
   if (~isempty(Tri1)),
      plot3(Tri1.XYZ(:,1), Tri1.XYZ(:, 2), Tri1.XYZ(:,3),'r*'); axis tight;hold on
   end
   plot3(ctr1(1), ctr1(2) , ctr1(3), 'b*');

   [err,PatchHandles] = ShowPlane (Eq1, Opt); view(3); hold on
   [err,PatchHandles] = ShowPlane (Eq2, Opt); view(3); hold on
   xlabel ('R'); ylabel ('A'); zlabel ('I');

   plot3(Line(:,1), Line(:,2), Line(:,3), 'r-', 'LineWidth', 4);
   stmp = sprintf('Rotation about intersection line from %s plane to oblique is %g degrees.\nMax displacement from cardinal plane %g mm\n', ...
                     CardPlane, AngleDeg, min(d));

   title (stmp);
   subplot 212; cla
   if (~isempty(Tri1)),
      plot3(Tri1.XYZ(:,1), Tri1.XYZ(:, 2), Tri1.XYZ(:,3),'r*'); axis tight;hold on
   end
   plot3(ctr1(1), ctr1(2) , ctr1(3), 'b*');

   [err,PatchHandles] = ShowPlane (Eq1, Opt); view(3); hold on
   [err,PatchHandles] = ShowPlane (Eq2, Opt); view(3); hold on
   xlabel ('Left to Right'); ylabel ('Posterior to Anterior'); zlabel ('Inferior to Superior');

   plot3(Line(:,1), Line(:,2), Line(:,3), 'r-', 'LineWidth', 4);
   axis equal; hold on


   %suggest a rotation with 3dvolreg, but first, make sure this is not a double oblique slice
   N_RotAx = 0
   if (Line(1,1) ~= Line (2,1)) ,
      RotAx(1) = 1;
      N_RotAx = N_RotAx + 1;
   else
      RotAx(1) = 0;
   end
   if (Line(1,2) ~= Line (2,2)) ,
      RotAx(2) = 1;
      N_RotAx = N_RotAx + 1;
   else
      RotAx(2) = 0;
   end
   if (Line(1,3) ~= Line (2,3)) ,
      RotAx(3) = 1;
      N_RotAx = N_RotAx + 1;
   else
      RotAx(3) = 0;
   end

   if (0 & N_RotAx > 1),
      fprintf (1,'Looks like you have super oblique slices, not ready to deal with that yet.\n');
   else
      if (RotAx(1)),
         AlignToOblique = sprintf ('3drotate -rotate %gR 0A 0I', AngleDeg);
      end
      if (RotAx(2)),
         AlignToOblique = sprintf ('3drotate -rotate 0R %gA 0I', AngleDeg);
      end
      if (RotAx(3)),
         AlignToOblique = sprintf ('3drotate -rotate 0R 0A %gI ', AngleDeg);
      end

      fidv(1) = 1;
      fidv(2) = fopen('DeOblique_out.txt','a');
      for (outp=1:1:2),
         %suggest 3drotate command
         fprintf (fidv(outp),'\nTo rotate Data to match %s:\n%s -prefix Data_Rotated_to_Oblique Data\nHit Enter to continue\n',...
                  fname, AlignToOblique);
         %if (outp == 1) pause; end
         fprintf (fidv(outp),'\nTypically, Data represents an anatomical \n');
         fprintf (fidv(outp),' data set acquired in a cardinal plane \n');
         fprintf (fidv(outp),' (Axial, Sagittal or Coronal). The oblique slice\n');
         fprintf (fidv(outp),' is from functional data. The suggested 3drotate\n');
         fprintf (fidv(outp),' command is meant to bring the anatomy (Data) into\n');
         fprintf (fidv(outp),' alignment with the function.\n');
         fprintf (fidv(outp),'If the rotation appears to be in the wrong direction,\n');
         fprintf (fidv(outp),' try using %g for a rotation angle.\n', -AngleDeg);
         fprintf (fidv(outp),' When that occurs, please notify me. \n\tsaadz@mail.nih.gov\n\n');
      end
   end

   end
   if (0), %other untested methods ....
      % ----------------------------------
      % get the rotation matrix
      %Making Obliques Cardinals:
      [err,ObliqueNowCard.XYZ, rot1] = AxisRotate3D(Tri1.XYZ, u , AngleRad, ctr1);
      for (i=1:1:3),
         fprintf (1,'Oblique p%g: | %.2f %.2f %.2f |\t%s p%g: | %.2f %.2f %.2f |\n', ...
               i, Tri1.XYZ(i,1), Tri1.XYZ(i,2), Tri1.XYZ(i,3), ...
               CardPlane, i, ObliqueNowCard.XYZ(i,1), ObliqueNowCard.XYZ(i,2), ObliqueNowCard.XYZ(i,3));
      end
      %Tri1.XYZ
      %ObliqueNowCard.XYZ
      %rot1

      %Making Cardinals Obliques
      [err,CardNowOblique.XYZ, rot2] = AxisRotate3D(ObliqueNowCard.XYZ, u , -AngleRad, ctr1);
      for (i=1:1:3),
         fprintf (1,'Cardinal p%g: | %.2f %.2f %.2f |\tObique p%g: | %.2f %.2f %.2f |\n', ...
                i, ObliqueNowCard.XYZ(i,1), ObliqueNowCard.XYZ(i,2), ObliqueNowCard.XYZ(i,3), ...
               i, CardNowOblique.XYZ(i,1), CardNowOblique.XYZ(i,2), CardNowOblique.XYZ(i,3));
      end

      %ObliqueNowCard.XYZ
      %CardNowOblique.XYZ
      %rot2

      %write results for 3drotate.
      % ************* SWAP rot1 and rot2 to make rtation go in the correct direction ...
      fmatvec = sprintf ('%s_To%s.matvec', fname, CardPlane);
      fid = fopen (fmatvec, 'w');
      if (fid < 0),
         fprintf (1,'Error: Failed to open output file %s for writing.\nCheck your permissions.\n', fmatvec);
      end
      for (i=1:1:3),
         fprintf (fid, '%g %g %g 0\n', rot2(i,1), rot2(i,2), rot2(i,3));
      end
      fclose (fid);
      fprintf (1,'To turn the oblique data to %s plane use:\n', CardPlane);
      fprintf (1,'3drotate -matvec_dicom %s ...\n\n',  fmatvec );


      [err, PathString, FileString] = GetPath(fname);
      fmatvec = sprintf ('%sCardinal_To%s.matvec', PathString, FileString);
      fid = fopen (fmatvec, 'w');
      if (fid < 0),
         fprintf (1,'Error: Failed to open output file %s for writing.\nCheck your permissions.\n', fmatvec);
      end
      for (i=1:1:3),
         fprintf (fid, '%g %g %g 0\n', rot1(i,1), rot1(i,2), rot1(i,3));
      end
      fclose (fid);
      fprintf (1,'To turn the cardinal data to %s''s plane use:\n', fname);
      fprintf (1,'3drotate -matvec_dicom %s ...\n\n',  fmatvec );

      fprintf (1,'To test for individual points add options: -origin %g %g %g -points\n\n',...
         ctr1(1), ctr1(2), ctr1(3));

   end

   fname = input ('Enter I filename (or test or dircos, or twovec): ', 's');
end
fprintf (fidv(2), '\n\n\tCreated with %s\n\t%s, %s\n', FuncName, pwd, date);
fclose(fidv(2));
