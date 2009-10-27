function app_neurofeedback(cfg)
% APP_NEUROFEEDBACK is a simple neurofeedback application which receives
% commands via an input stream and visualizes them using Psychtoolbox; it
% is designed to work with the bcifun_latidx and bcifun_frqclf examples
% functions. If the command is discrete then classification mode is
% assumed.
%
% Use as
%
%   app_neurofeedback(cfg)
%
% with the following configuration options:
%
% cfg.istream = the input stream that is used by read_event (default = []).
%               E.g., use 'tcp://localhost:1976' for TCP input on the presentation machine. 
% cfg.debug   = in debug mode we only give command line output (default = false)
%
% Copyright (C) 2009, Marcel van Gerven
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ~isfield(cfg,'istream'),  cfg.istream = []; end
if ~isfield(cfg,'debug'),    cfg.debug = false; end

% Make sure we run on OpenGL Psychtoolbox
if ~cfg.debug
    AssertOpenGL;
end

% number of read trials
n = 0;

% keep track of the sum
sumval = 0;

% keep track of the sum of squares
sumsqrval = 0;

try
  
  % Open onscreen window with default settings:
  if ~cfg.debug
  
      screenNumber = max(Screen('Screens'));
      [w,rect] = Screen('OpenWindow', screenNumber, [0 0 0], []);
      x = rect(3); y = rect(4); xc = x/2; yc = y/2; % center position
      
      HideCursor;
      
      % Draw fixation cross
      Screen('DrawLine', w, [100 100 100], xc-20, yc, xc+20, yc, 3);
      Screen('DrawLine', w, [100 100 100], xc, yc-20, xc, yc+20, 3);
      Screen('Flip',w,0,1);
      
  end
  
  while cfg.debug || ~KbCheck() 
  
    % get command via input stream
    event = read_event(cfg.istream);    
    
    if isfield(event,'value')
    
      cmd = event.value;
      fprintf('CMD: %s\n',num2str(cmd));
      
      if ~cfg.debug
        
        % compute colors
        
        if rem(cmd,1) % floats in neurofeedback mode
          
          % update mean and standard deviation
          n = n + 1;
          sumval = sumval + cmd;
          sumsqrval = sumsqrval + cmd.^2;
          
          % estimate the probability of deviation from the mean under a normal
          % distribution
          meanval = sumval ./ n;
          sdval = sqrt((sumsqrval - ((sumval.^2) ./ n)) ./ (n-1));
          
          p = 256 * max(min((cmd - meanval) ./ sdval,1),0);
          
        else % classification mode
          
          if cmd==1
            % left
            p = 0;
          else
            % right
            p = 256;
          end
          
        end
        
        Screen('FillOval', w, (256-p)*[1 1 1], [xc/2-50 yc-50 xc/2+50 yc+50]);
        Screen('FillOval', w, p*[1 1 1], [3*xc/2-50 yc-50 3*xc/2+50 yc+50]);
        
        Screen('Flip',w,0,1);
        
      end
    end
  end
  
catch

    if ~cfg.debug
        
        ShowCursor;
        Screen('CloseAll');
    
    end
    
    fprintf('%s\n',lasterr);

end

if ~cfg.debug
    ShowCursor;
    Screen('CloseAll');
end
