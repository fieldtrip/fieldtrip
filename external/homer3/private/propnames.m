function props = propnames(s)

if isstruct(s)
    props = fieldnames(s);
else
    props = properties(s);
end