%% 三个逻辑分布的DAEM模型函数
function F = daemo(x)
%x1-lnA x2-E1 x3-E2 x4-E3 x5-sigma1 x6-sigma2 x7-sigma3 x8-w1 x9-w2对应

%T0 = evalin('base', 'T0'); 
T = evalin('base', 'T'); 
R = 8.314e-3; % 气体常数,注意单位是 kJ/mol·K
beTa = 1/12; % 升温速率5 oC/min；
F = 0;
n = readtable('/Users/caoyang/Desktop/new_DAEM','ReadVariableNames',false);
%使用readTable读取数据文件
da = table2array(n(2:50,3))';

for i = 1:49 
f = @(E)(exp(x(1))./beTa).*exp(-E./(R.*T(i))-(exp(x(1))./beTa).*(R.*T(i).^2.*(...
    E+0.66691.*R.*T(i)))./(E.*(E+2.64943.*R.*T(i))).*exp(-E./(R.*T(i)))).*...
    (x(8)/(x(5).*(2*pi)^0.5).*exp(-(x(2)-E).^2/...
    (2*x(5)^2))+ x(9)/(x(6).*(2*pi)^0.5).*exp(-(x(3)-E).^2./...
    (2*x(6).^2))+ (1-x(8)-x(9))/(x(7).*(2*pi)^0.5).*exp(-(x(4)-E).^2./...
    (2*x(7).^2)));


dcal = quadgk(f,0,inf);
F = F+ (dcal - da(i))^2;
% da/dT 理论值和计算值之差平方和，累计到F上
end

end %函数到此结束