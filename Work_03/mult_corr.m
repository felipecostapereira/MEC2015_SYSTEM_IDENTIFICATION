function R2 = mult_corr(real,est)
% calculate R2 - multiple correlation coefficient

SSE = sum((real - est).^2);

avg_real = mean(real);

sum2 = sum((real - avg_real).^2);

R2 = 1 - SSE / sum2;

end