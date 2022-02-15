% Function check_value(var, vals, name)
%
% Checks the var (elements) to see if they belong to the values defined in
% vals
% returns false if it's in range true if not

function check_value(var, vals, name)
kk = find(ismember(var, vals) == 0);
if (~isempty(kk))
    error([ name ' may only contain the following values: ' num2str(vals) '.']);
end