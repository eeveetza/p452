function [emax, emin, mu_e_abs, mu_e, median_e, sigma_e, rmsd, rho] = predictionErrorStatistics(m,p)
% This function computes the error statistics of the predicted values as
% p(i) as compared to the measured ones e(i)
% Reference: 
% E. Ostlin, H. Suzuki, H-J. Zepernick, "Evaluation of the Propagation Model 
% Recommendation ITU-R P.1546 for Mobile Services in Rural Australia," 
%IEEE Trans on Vehicular Technology, Vol 57, No. 1, January 2008
%
% Rev   Date        Author                          Description
%-------------------------------------------------------------------------------
% v2    30AUG18     Ivica Stevanovic, OFCOM         Introduced median and RMSD
% v1    23AUG13     Ivica Stevanovic, OFCOM         Initial version

if (length(m) ~= length(p))
    warning('Vectors m and p must be of the same length')
    return
end

e=p-m;

emax=max((e));
emin=min((e));

mu_e_abs=mean(abs(e));

mu_e = mean(e);

sigma_e=std(e);

rho= sum((m-mean(m)).*(p-mean(p)))/(sqrt(sum((m-mean(m)).^2))*sqrt(sum((p-mean(p)).^2) ));

rmsd = sqrt((sum((m-p).^2))/length(m));

median_e = median(e);





