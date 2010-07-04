function afni_fig_interface(src, evnt)
   cf = gcf;
   figure(src);
   % evnt
   ud = get(src,'UserData');
   if (isempty(ud) | ~isfield(ud,'axes_lock')),
      ud = default_ud();
      set(src,'UserData',ud);
   end
   
   if (  lower(evnt.Character) == 'l'  ),
      ud.axes_lock = ~ud.axes_lock;
      set (src,'UserData',ud);
      if ud.axes_lock, set (src,'Name','L');
      else set (src,'Name','');
      end
   elseif (  lower(evnt.Character) == 'a'  ),
      if (ud.axes_lock),
         reset_axis(evnt.Character, get_all_axes(src));
      else 
         reset_axis(evnt.Character,gca); %can send arrays for characters and axes
      end 
   elseif (  lower(evnt.Character) == 'x'  |...
         lower(evnt.Character) == 'y'  |...
         lower(evnt.Character) == 'z' ),
      if (ud.axes_lock),
         zoom_axis(evnt.Character, get_all_axes(src));
      else 
         zoom_axis(evnt.Character,gca); %can send arrays for characters and axes
      end 
   elseif ( strcmp(evnt.Key , 'leftarrow') |...
            strcmp(evnt.Key , 'rightarrow') |...
            strcmp(evnt.Key , 'uparrow') |...
            strcmp(evnt.Key , 'downarrow')),
      if (ud.axes_lock),
         shift_axis(evnt, get_all_axes(src));
      else 
         shift_axis(evnt, gca);
      end 
   end
   figure(cf);
return;

function reset_axis(a,h)
   ca = gca;
   for (ih=1:1:length(h)),
      axes(h(ih));
      axis auto;
   end
   axes(ca);
return

function zoom_axis(a,h)
   ca = gca;
   for (ia=1:1:length(a)),
      for (ih=1:1:length(h)),
         axes(h(ih));
         dv = 4;
         if (  ~(a(ia) - lower(a(ia))) ),
            dv = 1;
         end
         if (~(lower(a(ia)) - 'x')), se = [1 2];
         elseif (~(lower(a(ia)) - 'y')), se = [3 4];
         elseif (~(lower(a(ia)) - 'z')), se = [5 6];
         else 
            return;
         end
         v = axis;
         if (length(v) < max(se)), return; end
         vc = (v(se(1))+v(se(2)))/2;
         vd = v(se(2))-v(se(1));
         v(se(1)) = vc-vd/dv;
         v(se(2)) = vc+vd/dv;
         axis(v);
      end   
   end
   axes(gca);
return;

function shift_axis(e,h)
   ca = gca;
   for (ie=1:1:length(e)),
      for (ih=1:1:length(h)),
         axes(h(ih));
         if (strcmp(e(ie).Key , 'leftarrow')), 
            se = [1 2]; dv = +0.1;
         elseif (strcmp(e(ie).Key , 'rightarrow')), 
            se = [1 2]; dv = -0.1;
         elseif (strcmp(e(ie).Key , 'downarrow')), 
            se = [3 4]; dv = -0.1;
         elseif (strcmp(e(ie).Key , 'uparrow')), 
            se = [3 4]; dv = 0.1;
         else 
            return;
         end
         v = axis
         vd = v(se(2))-v(se(1));
         v(se(1)) = v(se(1))+vd*dv;
         v(se(2)) = v(se(2))+vd*dv;
         axis(v);
      end   
   end
   axes(gca);
return;

function v = get_all_axes(f)
   l = get(f, 'Children');
   v = [];
   for (i=1:1:length(l)),
      if (strcmp(get(l(i),'Type'),'axes')), v = [v l(i)]; end
   end
return

function ud = default_ud()
   ud = struct('axes_lock', 0);
return
