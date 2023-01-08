%%新函数设计，假设指前因A相同
function F = Pseudo_DAEM3(x)
R = 8.314e-3; % Gas constant 8.314*10-3 kJ/mol*K;
beTa = 1/6; % Heating rate, 10 K/min;
da =  load("/Users/caoyang/Desktop/DAEM/20wt-0.2/da.txt"); 
T = evalin('base', 'T'); 
da = da'; 
SSE = 0;
SST = 0;
mea = mean(da);
da = da(1:10:end);
% Your data
% Transpose of da
% x(1) log_10^A x(2) E1 x(3) sigma1 x(4) E2 x(5) sigma2 x(6) E3 x(7) sigma3 
% x(8) E4 x(9) sigma4 x(10) w1 x(11)w2 x(12) w3 
% Pre-expoential factor ranged from 10^10 to 10^20
for i = 1:285
fun = @(E)...
10.^17.5/beTa.*(exp(...
    -E./(R.*T(i))-(10.^17.5./beTa).*(T(i).*exp(-E./(R.*T(i)))-E./R.*expint(E./(R.*T(i))) ...
    ))).*(x(7).*...
(1./((2.*pi).^0.5.*x(2))).*exp(-(E-x(1)).^2/(2.*x(2).^2))...
+x(8).*...
(1./((2.*pi).^0.5.*x(4))).*exp(-(E-x(3)).^2/(2.*x(4).^2)) ...
+x(9).*...
(1./((2.*pi).^0.5.*x(6))).*exp(-(E-x(5)).^2/(2.*x(6).^2)) ...
);
dcal =  quadgk(fun,0,500); %使用MATLAB的自适应积分函数计算
SSE = SSE+ (dcal - da(i))^2;
SST = SST + (mea-da(i))^2;
end
F= SSE./SST;
end
