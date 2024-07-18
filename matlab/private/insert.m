function y = insert(x, idx, a)
%insert(x, idx, a) inserts elements from vector a at a position idx of vector x

% Rev   Date        Author                     Description
%-------------------------------------------------------------------------------
% v1    18JUL23     Ivica Stevanovic, OFCOM    Initial version


% make sure x and a are row vectors
[nr, nc] = size (x);
if min(nr, nc) >  1
    error('x must be a row vector')
end

if(nr > 1)
    x = x.';
end

[nr, nc] = size (a);
if min(nr, nc) >  1
    error('a must be a row vector or a scalar')
end

if(nr > 1)
    a = a.';
end

y = [x(1:end<idx), a, x(1:end >= idx)];

return
end
