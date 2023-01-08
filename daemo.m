%% 数f据读取及常量设置
clear;clc;
func = @Pseudo_DAEM3;
T = load("/Users/Desktop/DAEM/20wt-0.2/T.txt");
%T所在的txt文件，同函数，尽量导出txt数据（制表符分隔）不要直接用xls（用超算请忽略）
T = T';
T = T(1:10:end);
%% 模式搜索部分
x0 = [150, 15, 190, 5, 220, 35, 0.3, 0.45, 0.25];
Aeq = [0 0 0 0 0 0 1 1 1];
beq = 1; % w1+w2+w3=1 
A = [];     
b = [];
lb = [100, 1.5, 100, 1.5, 100, 1.5, 0, 0, 0]; %变量的下限
ub = [350, 75, 350, 75, 350, 75, 1, 1, 1]; %变量的上限
nonlcon = [];   
options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf,'MaxIter',...
15000, 'TolMesh',5e-4);
% 条件设置：模式搜索法，迭代次数，最大迭代次数15000， TolMesh为网格尺寸阈值。
x = patternsearch(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%% 三个伪组分的Dalpha/DT的热解动力学模型；
R = 8.314e-3; % 气体常数8.314*10-3 kJ/mol*K;
beTa = 1/6; % 升温速率10 K/min;
alp = [x(1); x(2); x(7)];  % Pseudo component1 E1, sigma1, w1;
bet = [x(3); x(4); x(8)];  % Pseudo component2 E2, sigma2, w2;
gam = [x(5); x(6); x(9)];  % Pseudo component3 E3, sigma3, w3;
for i = 1:285 
d1 = @(E) 10.^17.5/beTa.*(exp(...
    -E./(R.*T(i))-(10.^17.5./beTa).*(T(i).*exp(-E./(R.*T(i)))-E./R.*expint(E./(R.*T(i)))))).*alp(3).*(1./((2.*pi).^0.5.*alp(2))).*exp(-(E-alp(1)).^2/(2.*alp(2).^2));
d2 = @(E) 10.^17.5/beTa.*(exp(...
    -E./(R.*T(i))-(10.^17.5./beTa).*(T(i).*exp(-E./(R.*T(i)))-E./R.*expint(E./(R.*T(i)))))).*bet(3).*(1./((2.*pi).^0.5.*bet(2))).*exp(-(E-bet(1)).^2/(2.*bet(2).^2));
d3 = @(E) 10.^17.5/beTa.*(exp(...
    -E./(R.*T(i))-(10.^17.5./beTa).*(T(i).*exp(-E./(R.*T(i)))-E./R.*expint(E./(R.*T(i)))))).*gam(3).*(1./((2.*pi).^0.5.*gam(2))).*exp(-(E-gam(1)).^2/(2.*gam(2).^2));
D1(i) = quadgk(d1,0,500);
D2(i) = quadgk(d2,0,500);
D3(i) = quadgk(d3,0,500);
% D1~D3为伪组分的DTG绘图数据
end  
%% 总拟合曲线及伪组分拟合数据生成
R_esult = x
T = T';
func2= @da_dT3;
dadT = func2(x);
dadT = dadT';
D1 = D1';
D2 = D2';
D3 = D3';
