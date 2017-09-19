%% 
function Ct = dualinput_nlls(coeff,xdata)

a_coeff = coeff(1);
p_coeff = coeff(2);
k_coeff = coeff(3);

aorta = xdata(:,1);
portal = xdata(:,2);
T = xdata(:,3);

n = numel(T);
Ct = zeros(n,1);
dt=T(2)-T(1);
F0 = (a_coeff*aorta+p_coeff*portal)';
F1 = repmat(F0, [n 1]).*tril(toeplitz(exp(-k_coeff*(T-T(1)))));
Ct = (T(2)-T(1))*(sum(F1,2) - 0.5*(F1(:,1)+diag(F1)));
Ct(1)=0;
