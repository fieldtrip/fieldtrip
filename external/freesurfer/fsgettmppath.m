function tmppath = fsgettmppath(tmppathdefault)
% tmppath = fsgettmppath(<tmppathdefault>)
% Gets path to a temporary folder (does NOT generate a random file name)
% First looks in the following order:
%  1. $TMPDIR - env var must exist and folder must exist
%  2. $TEMPDIR - env var must exist and folder must exist
%  3. /tmp - folder must exist   -> JMS changed order here
%  4. /scratch - folder must exist -> JMS changed order here
%  5. tmppathdefault - if passed
%  6. current folder, ie, ./ (prints warning)


tmppath = getenv('TMPDIR');
if(~isempty(tmppath)) 
  if(exist(tmppath)) 
    return; 
  end
end

tmppath = getenv('TEMPDIR');
if(~isempty(tmppath)) 
  if(exist(tmppath)) 
    return; 
  end
end  

% JMS: let /tmp prevail, rather than scratch, at DCCN we don't have write
% permission in /scratch
tmppath = '/tmp';
if(exist(tmppath)) return; end

tmppath = '/scratch';
if(exist(tmppath)) return; end


if(nargin > 0)
  tmppath = tmppathdefault;
  if(exist(tmppath)) return; end
end

tmppath = './';
fprintf(['WARNING: fsgettmppath: could not find a temporary folder,' ...
	 ' using current folder\n']);

return;



