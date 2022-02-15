% Function check_limit(var, low, hi, name, type)
%
% Checks the var to see if it's inbetween the range low <= var <= hi
% returns false if it's in range true if not

function check_limit(var, low, hi, name)
if ((var < low) || (var > hi))
    error([ name ' = ' num2str(var) ' is outside the limits: [' num2str(low) ', ' num2str(hi) '].']);
end
