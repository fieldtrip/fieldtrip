function IntProj(ana, hs, proj)
%
%  A function to perform Minimum or Maximum intensity projection of a volume of data
%  ana: Name of AFNI brick file
%  hs: Number of slices to use is each direction. Projection at slice i is from i-hs to i+hs
%      hs stands for half slab
%  proj: (-1) minimum intensity projection
%          1  maximum intensity projection
%  For output, the function creates 
% see also insane function Int_Proj

FuncName = 'IntProj';

if (proj == -1) proj_s = 'min';
elseif (proj == 1) proj_s = 'max';
else 
   fprintf(2,'Error %s: Bad projection parameter\n', FuncName);
   return;
end

[ana_pref, ana_view] = AfniPrefix(ana);
ana = sprintf('%s%s', ana_pref, ana_view);

[err, V, V_Info] = BrikLoad(ana);
Nvox = V_Info.DATASET_DIMENSIONS;

verb = 1;

%padd V by half slab in all directions to fix edge effects
Vp = zeros(Nvox(1)+3.*hs, Nvox(2)+3.*hs, Nvox(3)+3.*hs);
Vp(hs+1:Nvox(1)+hs, hs+1:Nvox(2)+hs, hs+1:Nvox(3)+hs) = V; clear V;


for (pd = 1:1:3), % projection direction
   if (pd == 1) dc = 'I';
   elseif (pd == 2) dc = 'J';
   elseif (pd == 3) dc = 'K';
   else 
      fprintf(2,'Error %s: Bad projection direction\n', FuncName);
      return;
   end
   
   if (verb) fprintf(1,'%s: Processing direction %s\n', FuncName, dc); end
   %create the output set
   Vmi = zeros(Nvox(1), Nvox(2), Nvox(3));
   
   if (proj == -1),
      if (pd == 1),
         for (i=hs+1:1:Nvox(1)+hs),
            Vmi(i-hs,:,:) = min(Vp(i-hs:i+hs,hs+1:Nvox(2)+hs,hs+1:Nvox(3)+hs), [], pd);
         end
         Vmall = Vmi;	
      elseif (pd == 2),
         for (i=hs+1:1:Nvox(2)+hs),
            Vmi(:,i-hs,:) = min(Vp(hs+1:Nvox(1)+hs,i-hs:i+hs,hs+1:Nvox(3)+hs), [], pd);
         end
         Vmall(:) = min(Vmi(:),Vmall(:));	
      elseif (pd == 3),
         for (i=hs+1:1:Nvox(3)+hs),
            Vmi(:,:,i-hs) = min(Vp(hs+1:Nvox(1)+hs,hs+1:Nvox(2)+hs,i-hs:i+hs), [], pd);
         end	
         Vmall(:) = min(Vmi(:),Vmall(:));	
      end
   elseif (proj == 1),
      if (pd == 1),
         for (i=hs+1:1:Nvox(1)+hs),
            Vmi(i-hs,:,:) = max(Vp(i-hs:i+hs,hs+1:Nvox(2)+hs,hs+1:Nvox(3)+hs), [], pd);
         end
         Vmall = Vmi;	
      elseif (pd == 2),
         for (i=hs+1:1:Nvox(2)+hs),
            Vmi(:,i-hs,:) = max(Vp(hs+1:Nvox(1)+hs,i-hs:i+hs,hs+1:Nvox(3)+hs), [], pd);
         end
         Vmall(:) = max(Vmi(:),Vmall(:));	
      elseif (pd == 3),
         for (i=hs+1:1:Nvox(3)+hs),
            Vmi(:,:,i-hs) = max(Vp(hs+1:Nvox(1)+hs,hs+1:Nvox(2)+hs,i-hs:i+hs), [], pd);
         end	
         Vmall(:) = max(Vmi(:),Vmall(:));	
      end
   
   end
   OptW.Prefix = sprintf('%s_%s%g%s', ana_pref, proj_s, hs, dc);
   OptW.Scale = 1;
   WriteBrik(Vmi, V_Info, OptW);
end

   OptW.Prefix = sprintf('%s_%s%gIJK', ana_pref, proj_s, hs);
   OptW.Scale = 1;
   WriteBrik(Vmall, V_Info, OptW);

return;
