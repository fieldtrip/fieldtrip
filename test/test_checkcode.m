function test_checkcode

% WALLTIME 00:20:00
% MEM 2gb

% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3309

%%

[v, p] = ft_version;
p = {p};
f = {};

while ~isempty(p)
  fprintf('looking in directory %s\n', p{1});
  nf = dir(p{1});
  p(1) = []; % remove this one
  for i=1:numel(nf)
    if nf(i).isdir && ~isequal(nf(i).name(1), '.') && ~isequal(nf(i).name, 'external') && ~isequal(nf(i).name, 'test')
      p{end+1} = fullfile(nf(i).folder, nf(i).name);
    else
      if ~nf(i).isdir && nf(i).name(end-1)=='.' && nf(i).name(end)=='m'
        f{end+1} = fullfile(nf(i).folder, nf(i).name);
      end
    end
  end
end

filelist = f;
clear f p

fprintf('found %d *.m files\n', numel(filelist));

%%

invalid = {
  'Use || instead of | as the OR operator in (scalar) conditional statements.'
  'Use && instead of & as the AND operator in (scalar) conditional statements.'
  'is known to MATLAB by its file name'
  'Use ISCELL instead of comparing the class to ''cell''.'
  'Use ISLOGICAL instead of comparing the class to ''logical''.'
  'Use TRUE or FALSE instead of LOGICAL(1) or LOGICAL(0).'
  };

m = {};
for i=1:numel(filelist)
  f = filelist{i};
  s = checkcode(f);
  m = union(m, {s.message});
  for j=1:numel(s)
    for k=1:numel(invalid)
      if contains(s(j).message, invalid{k})
        % pretty display
        fprintf('================================================================================\n');
        fprintf('ERROR: %s\n', invalid{k});
        fprintf('================================================================================\n');
        
        checkcode(f)
        error('please fix line %d in %s', s(j).line, f);
      end
    end
  end
end

%%

% this shows the full list of warning messages
if false
  disp(m);
end

% this is a subset of the LINT warnings, I removed the ones about unused variables
%
%     ''findstr' is not recommended. Use 'strfind' instead.'
%     ''flipdim' is not recommended. Use 'flip' instead.'
%     ''isempty(strfind(str1, str2))' is not recommended. Use '~contains(str1, str2)' instead.'
%     ''isequalwithequalnans' is not recommended. Use 'isequaln' instead.'
%     ''resampled_tsdif' produces a value that might be unused.'
%     ''skip' produces a value that might be unused.'
%     ''spectrum.welch' is not recommended. Use 'pwelch' instead.'
%     ''strfind(str1, str2)' is not recommended. Use 'contains(str1, str2)' instead.'
%     ''strread' is not recommended. Use 'textscan' instead.'
%     ''textread' is not recommended. Use 'textscan' instead.'
%     ''||' produces a value that might be unused.'
%     ''~cellfun('isempty', strfind(str1, str2))' is not recommended. Use 'contains(str1, str2)' instead.'
%     ''~cellfun(@isempty, strfind(str1, str2))' is not recommended. Use 'contains(str1, str2)' instead.'
%     ''~isempty(strfind(str1, str2))' is not recommended. Use 'contains(str1, str2)' instead.'
%     'Assign the onCleanup output argument to a variable. Do not use the tilde operator (~) in place of a variable.'
%     'Best practice is to separate output variables with commas.'
%     'Calling AXES(h) in a loop can be slow. Consider moving the call to AXES outside the loop.'
%     'Consider replacing GRIDDATA with SCATTEREDINTERPOLANT for better performance.'
%     'Consider using ISA instead of comparing the class name.'
%     'Consider using newline, semicolon, or comma before this statement for readability.'
%     'DISP(SPRINTF(...)) can usually be replaced by FPRINTF(...\n).'
%     'ERROR takes SPRINTF-like arguments directly.'
%     'EXIST with two input arguments is generally faster and clearer than with one input argument.'
%     'Extra comma is unnecessary in CASE statement before newline.'
%     'Extra comma is unnecessary in CATCH statement before newline.'
%     'Extra comma is unnecessary in ELSE statement before newline.'
%     'Extra comma is unnecessary in ELSEIF statement before newline.'
%     'Extra comma is unnecessary in END statement before newline.'
%     'Extra comma is unnecessary in FOR statement before newline.'
%     'Extra comma is unnecessary in IF statement before newline.'
%     'Extra comma is unnecessary in OTHERWISE statement before newline.'
%     'Extra comma is unnecessary in SWITCH statement before newline.'
%     'Extra comma is unnecessary in WHILE statement before newline.'
%     'Extra comma is unnecessary.'
%     'Extra semicolon is unnecessary in CASE statement before newline.'
%     'Extra semicolon is unnecessary in ELSE statement before newline.'
%     'Extra semicolon is unnecessary in ELSEIF statement before newline.'
%     'Extra semicolon is unnecessary in END statement before newline.'
%     'Extra semicolon is unnecessary in FOR statement before newline.'
%     'Extra semicolon is unnecessary in FUNCTION statement before newline.'
%     'Extra semicolon is unnecessary in IF statement before newline.'
%     'Extra semicolon is unnecessary in WHILE statement before newline.'
%     'Extra semicolon is unnecessary.'
%     'FOR might not be aligned with its matching END (line 169).'
%     'FOR might not be aligned with its matching END (line 31).'
%     'FOR might not be aligned with its matching END (line 437).'
%     'FOR might not be aligned with its matching END (line 81).'
%     'FREAD(FID,...,'*char') is more efficient than CHAR(FREAD(...)).'
%     'For improved robustness, consider replacing i and j by 1i.'
%     'For readability and performance, consider using 'newline' instead of 'sprintf('\n')'.'
%     'Function ERROR might be called with too few arguments.'
%     'IF might not be aligned with its matching END (line 172).'
%     'IF might not be aligned with its matching END (line 30).'
%     'IF might not be aligned with its matching END (line 37).'
%     'IF might not be aligned with its matching END (line 54).'
%     'INV(A)*b can be slower and less accurate than A\b. Consider using A\b for INV(A)*b or b/A for b*INV(A).'
%     'Instead of using transpose ('), consider using a different DIMENSION input argument to MAX.'
%     'Instead of using transpose ('), consider using a different DIMENSION input argument to MEAN.'
%     'Instead of using transpose ('), consider using a different DIMENSION input argument to MIN.'
%     'Instead of using transpose ('), consider using a different DIMENSION input argument to SUM.'
%     'LASTERR and LASTERROR are not recommended. Use an identifier on the CATCH block instead.'
%     'Loop index 'i' is changed inside of a FOR loop.'
%     'NARGCHK is not recommended. Use NARGINCHK without ERROR instead.'
%     'NUMEL(x) is usually faster than PROD(SIZE(x)).'
%     'Operator '+' is seldom used in a logical context.'
%     'Overloading DISPLAY is not recommended.'
%     'Parentheses are not needed in a FOR statement.'
%     'Parse error at FOR: usage might be invalid MATLAB syntax.'
%     'Possible inappropriate use of - operator. Use = if assignment is intended.'
%     'RAND or RANDN with the 'seed', 'state', or 'twister' inputs is not recommended. Use RNG instead.'
%     'STRMATCH is not recommended. Use STRCMP instead.'
%     'STRMATCH is not recommended. Use STRNCMP or VALIDATESTRING instead.'
%     'TRY might not be aligned with its matching END (line 299).'
%     'TRY statement should have a CATCH statement to check for unexpected errors.'
%     'Terminate statement with semicolon to suppress output (in functions).'
%     'Terminate statement with semicolon to suppress output (within a script).'
%     'The comparison will likely fail due to case mismatch.'
%     'The format might not agree with the argument count.'
%     'The warning with tag 'MATLAB:divideByZero' has been removed from MATLAB, so this statement has no effect.'
%     'This sparse indexing expression is likely to be slow.'
%     'This statement (and possibly following ones) cannot be reached.'
%     'This use of MAT2CELL should probably be replaced by a simpler, faster call to NUM2CELL.'
%     'To improve performance, replace ISEMPTY(FIND(X)) with ISEMPTY(FIND( X, 1 )).'
%     'Unexpected use of ' in a scalar context.'
%     'Unexpected use of '.*' in a scalar context.'
%     'Unexpected use of '.^' in a scalar context.'
%     'Unexpected use of '[' in a scalar context.'
%     'Use FIND with the 'first' or 'last' option.'
%     'Use STRCMP instead of == or ~= to compare character vectors.'
%     'Use STRCMPI(str1,str2) instead of using UPPER/LOWER in a call to STRCMP.'
%     'Use STRTRIM(str) instead of nesting FLIPLR and DEBLANK calls.'
%     'Use dynamic fieldnames with structures instead of GETFIELD.'
%     'Use dynamic fieldnames with structures instead of SETFIELD.'
%     'Use of 'nargout' in a script will be removed in a future release.'
%     'Use of brackets [] is unnecessary. Use parentheses to group, if needed.'
%     'Use one call to TEXTSCAN instead of calling STRTOK in a loop.'
%     'Using ISEMPTY is usually faster than comparing LENGTH to 0.'
%     'Using REGEXP(str, pattern, 'ONCE') is faster in this case.'
%     'Variable 'ElecPres' might be set by a nonlogical operator.'
%     'Variable 'ElecPres' might be set by a nonscalar operator.'
%     'Variable 'H', apparently a structure, is changed but the value seems to be unused.'
%     'Variable 'Npoints' might be set by a nonscalar operator.'
%     'Variable 'Nyq' might be set by a nonscalar operator.'
%     'Variable 'binwidth' might be set by a nonscalar operator.'
%     'Variable 'cfg', apparently a structure, is changed but the value seems to be unused.'
%     'Variable 'cond' might be set by a nonlogical operator.'
%     'Variable 'cond' might be set by a nonscalar operator.'
%     'Variable 'hdr', apparently a structure, is changed but the value seems to be unused.'
%     'Variable 'info', apparently a structure, is changed but the value seems to be unused.'
%     'Variable 'locktrllen' might be set by a nonscalar operator.'
%     'Variable 'nchans' might be set by a nonscalar operator.'
%     'Variable 'ncmb' might be set by a nonscalar operator.'
%     'Variable 'npad' might be set by a nonscalar operator.'
%     'Variable 'origrank' might be set by a nonscalar operator.'
%     'Variable 'surface', apparently a structure, is changed but the value seems to be unused.'
%     'WARNING takes SPRINTF-like arguments directly.'
%     'When checking if a variable is a matrix consider using ISMATRIX.'
%     '{ A{:} B } can often be replaced by [ A {B}], which can be much faster.'


