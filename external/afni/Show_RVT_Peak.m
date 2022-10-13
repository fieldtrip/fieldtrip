function Show_RVT_Peak(R,fg),
   for (icol = 1:1:length(R)),
      figure(fg); clf
      set (fg, 'KeyPressFcn', @afni_fig_interface);
      subplot(211);
      plot (R(icol).t, real(R(icol).X),'g'); hold on
      if (~isempty(R(icol).RVT)),
         subplot (211);
            plot (R(icol).tmidprd,...
                  zscale(R(icol).RVT, max(R(icol).ptrace),...
                  min(R(icol).ptrace)), 'k');
      end
      plot( R(icol).tptrace, R(icol).ptrace,'ro',...
            R(icol).tptrace, R(icol).ptrace,'r');
      plot( R(icol).tntrace, R(icol).ntrace,'bo',...
            R(icol).tntrace, R(icol).ntrace,'b');
      plot (R(icol).tmidprd, R(icol).ptracemidprd,'kx');
      for (i=1:1:length(R(icol).prd)),
       text( R(icol).tmidprd(i), R(icol).ptracemidprd(i),...
             sprintf('%.2f', R(icol).prd(i)));
      end
      if (isfield(R(icol), 'tR')),
         if (isfield(R(icol), 'ptraceR')),
            plot( R(icol).tR, R(icol).ptraceR,'m');
            plot( R(icol).tR, R(icol).ntraceR,'y');
         end
         if (  isfield(R(icol), 'RVTRS')),
               plot( R(icol).tR,...
                     zscale(  R(icol).RVTRS, max(R(icol).ptrace),...
                              min(R(icol).ptrace) ),...
                     'k.');
         end
      end
      xlabel('time (sec)');
      title (R(icol).vname, 'Interpreter', 'None');

      subplot (413);
      vn = real(R(icol).X)./(abs(R(icol).X)+eps);
      plot (R(icol).t, vn, 'g'); hold on

      plot (R(icol).t, R(icol).phz./2./pi, 'm');
      if (isfield(R(icol), 'phzR')),
         plot (R(icol).tR, R(icol).phzR./2./pi, 'm-.');
      end
      xlabel('time (sec)');
      title ('Scaled by magnitude of analytical signal', 'Interpreter', 'None');
      legend({'Scaled signal','phase'});

      subplot (414);
      plot (R(icol).tst, R(icol).phz_slc(:,1), 'ro', ...
            R(icol).tst, R(icol).phz_slc(:,2), 'bo', ...
            R(icol).tst, R(icol).phz_slc(:,2), 'b-'); hold on;
      plot (R(icol).t, R(icol).phz, 'k');
      grid on
      xlabel('time (sec)');
      title ('Phase sampled at slice acquisition time');
      legend({'slice 0', 'slice 1', 'slice 1', 'original phase'});
      plotsign2(fg);
      drawnow;
   end
return;
