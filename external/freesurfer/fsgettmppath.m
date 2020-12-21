function tmppath = fsgettmppath(tmppathdefault)
% tmppath = fsgettmppath(<tmppathdefault>)
% Gets path to a temporary folder (does NOT generate a random file name)
% First looks in the following order:
%  1. $TMPDIR - env var must exist and folder must exist
%  2. $TEMPDIR - env var must exist and folder must exist
%  3. /scratch - folder must exist
%  4. /tmp - folder must exist
%  5. tmppathdefault - if passed
%  6. current folder, ie, ./ (prints warning)


tmppath = getenv('TMPDIR');
if(~isempty(tmppath)) 
  if(exist(tmppath) == 7) 
    return; 
  end
end

tmppath = getenv('TEMPDIR');
if(~isempty(tmppath)) 
  if(exist(tmppath) == 7) 
    return; 
  end
end  

tmppath = '/scratch';
if(exist(tmppath) == 7) return; end

tmppath = '/tmp';
if(exist(tmppath) == 7) return; end

if(nargin > 0)
  tmppath = tmppathdefault;
  if(exist(tmppath) == 7) return; end
end

tmppath = './';
fprintf(['WARNING: fsgettmppath: could not find a temporary folder,' ...
	 ' using current folder\n']);

return;
