function b = isproperty(obj, propname)

if isstruct(obj)
    b = isfield(obj, propname);
elseif isobject(obj)
    b = all(isprop(obj, propname));
else
    b = false;
end
