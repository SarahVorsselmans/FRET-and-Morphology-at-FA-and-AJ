function tf = isStringScalar(x)
    tf = isstring(x) && isscalar(x);
end