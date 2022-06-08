%% 模式搜索部分(适用于3个密度函数分布)
clear;clc;
fun = @daemo;
%daTa = load('')
%T0 = daTa(1);
%T = daTa(2:40);
n = readtable('/Users/caoyang/Desktop/new_DAEM','ReadVariableNames',false);
%n为读取数据的表格
Tx  = table2array(n(1:50,1));
T0 = Tx(1);
T = Tx(2:50)';

options=optimoptions('patternsearch','Display','iter','PlotFcn',...
    @psplotbestf,'MaxIter',50000,...
    'MaxFunEvals',200000,'TolMesh',1e-7,'TolX',1e-15);
% 模式搜索法，迭代次数，最大迭代次数50000

%% 变量以及初值设置
x0 = [31 152 183 203 3.8 5.4 9.6 0.24 0.35]; %初始点随便选啦，只要在范围内
% x = [lnA E1 E2 E3 sigma1 sigma2 sigma3 w1 w2]
A = [0 0 0 0 0 0 0 1 1]; % 妹说就是零
b = [1];
Aeq = [];
beq =[];
lb = [20 75 75 75 0.5 0.5 0.5 0 0 ];% 变量的下限
ub =[45 300 300 300 50 50 50 1 1];% 变量的上限
nonlcon = [];
x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);


