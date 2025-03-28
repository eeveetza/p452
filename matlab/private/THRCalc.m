function thr = THRCalc(m,p,lt)
% This function computes total hit rate of prediction p as compared to the
% measurements m at different threshold values lt
%
% Rev   Date        Author                          Description
%-------------------------------------------------------------------------------
% v1    23AUG13     Ivica Stevanovic, OFCOM         Initial version

if (length(m) ~= length(p))
    warning('Vectors m and p must be of the same length')
    return
end

N=length(m);

for ii=1:length(lt)
    thr(ii)=( sum( (m>=lt(ii)).*(p>=lt(ii)) ) + ...
              sum( (m <lt(ii)).*(p <lt(ii))) )/N;
end
end

