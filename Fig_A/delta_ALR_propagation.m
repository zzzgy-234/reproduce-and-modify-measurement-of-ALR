clear;
% 计算亮度校准中，P+与P-的关系
syms A_LR P_plus P_minu s1 s2 s3 s4 e d e1 d1 delt_e delt_d delta_e1 N1 N2 N3 N4%定义参数

% 定义函数 f_alr
f_alr = @(s1, s2, s3, s4) (4 .* (s1 - s3) .* (s2 - s3)) ./ (s3 .* (s1 - s2)) - (4 .* (s1 - s3) .* (s2 - s3) .* ((s1 .* (-s3 ./ (4 .* (s1 - s3) .* (s2 - s3) .* (s1 + s2 - s3 - s4))).^(1/2)) ./ 2 - (s2 .* (-s3 ./ (4 .* (s1 - s3) .* (s2 - s3) .* (s1 + s2 - s3 - s4))).^(1/2)) ./ 2 + 1)) ./ (s3 .* (s1 - s2));

% 计算f_alr的偏导数并平方和
grad_f_alr = gradient(f_alr(s1, s2, s3, s4), [s1, s2, s3, s4]);
disp(grad_f_alr);
theo_delta_alr = sqrt((grad_f_alr(1).*sqrt(1./N1).*s3.*(1-P_minu.*A_LR)) .^ 2+(grad_f_alr(2).*sqrt(1./N2).*s3.*(1+P_plus.*A_LR)) .^ 2+(grad_f_alr(3).*sqrt(1./N3).*s3) .^ 2+(grad_f_alr(4).*sqrt(1./N4).*s3.*(1-P_plus*P_minu+(P_plus-P_minu).*A_LR)) .^ 2);

%进行符号替换
theo_delta_alr = subs(theo_delta_alr, {'s1','s2','s4'}, {str2sym('s3*(1-P_minu*A_LR)'), str2sym('s3*(1+P_plus*A_LR)'), str2sym('s3*(1-P_plus*P_minu+(P_plus-P_minu)*A_LR)')});

disp(theo_delta_alr);
% 定义数值
P_plus_val = linspace(0.1,0.8,1000);
P_minu_val = linspace(0.1, 0.8, 1000);
s3_val = 4.0;
A_LR_val = 0.16;
N1_val=0.25*3*10^6;
N2_val=0.25*3*10^6;
N3_val=0.25*3*10^6;
N4_val=0.25*3*10^6;

% 将符号表达式转换为数值函数
theo_erro_alr_numeric = matlabFunction(theo_delta_alr, 'Vars', {P_minu,P_plus, s3, A_LR, N1, N2,N3,N4});

% 计算数值表达式 
theo_delta_alr_values1 = theo_erro_alr_numeric(P_minu_val, P_plus_val,s3_val,A_LR_val, N1_val, N2_val,N3_val,N4_val);


% 定义符号变量
syms s1 s2 P_minu A_LR S_u N1 N2 dp

% 标准方案
eqns = [
    s1 == S_u * (1 - P_minu * A_LR);
    s2 == S_u * (1 + P_minu * A_LR);
];

[A_LR_sol, S_u_sol] = solve(eqns, [A_LR, S_u]);

% 显示结果
disp('A_LR = ');
disp(A_LR_sol);
disp('S_u=');
disp(S_u_sol);

% 定义函数 f_alr
f_alr = @(s1, s2, P_minu) -(s1 - s2)./(P_minu.*(s1 + s2));

% 计算理论误差
grad_f_alr = gradient(f_alr(s1, s2, P_minu), [s1, s2, P_minu]);
theo_erro_alr = sqrt((grad_f_alr(1).*sqrt(1./N1).*S_u.*(1-P_minu.*A_LR)).^2 + (grad_f_alr(2).*sqrt(1./N2).*S_u.*(1+P_minu.*A_LR)).^2 + (grad_f_alr(3).*(dp).*P_minu).^2);

% 进行符号替换
theo_erro_alr = subs(theo_erro_alr, {'s1','s2'}, {S_u.*(1-P_minu.*A_LR), S_u.*(1+P_minu.*A_LR)});

% 显示新的表达式
disp('the theoretical alr error =');
disp(theo_erro_alr);

% 变量赋值
A_LR_val = 0.16;
N1_val = 0.5 * 3 * 10^6;
N2_val = 0.5 * 3 * 10^6;
P_minu_val = linspace(0.1, 0.8, 1000);
S_u_val = 2; % 将符号值转换为双精度值

% 将符号表达式转换为数值函数
theo_erro_alr_numeric = matlabFunction(theo_erro_alr, 'Vars', {P_minu, S_u, A_LR, dp, N1, N2});

% 计算数值表达式 dp = 0.01
dp_val1 = 0.01;
theo_erro_alr_values1 = theo_erro_alr_numeric(P_minu_val, S_u_val, A_LR_val, dp_val1, N1_val, N2_val);

% 计算数值表达式 dp = 0.02
dp_val2 = 0.02;
theo_erro_alr_values2 = theo_erro_alr_numeric(P_minu_val, S_u_val, A_LR_val, dp_val2, N1_val, N2_val);

% 绘图
figure;
plot(P_minu_val, theo_delta_alr_values1, 'g-');
hold on;
plot(P_minu_val, theo_erro_alr_values1, 'b-', 'DisplayName', 'dp = 0.01');
hold on;
plot(P_minu_val, theo_erro_alr_values2, 'r--', 'DisplayName', 'dp = 0.02');
hold off;

% 添加文本注释
txt1 = {'\deltaP/P=1%'};
txt2 = {'\deltaP/P=2%'};
txt3 = {'A_{LR}=0.16'};
txt4={'3\bullet10^{6} Z evets'};
txt5={'4 bunch scheme'};
text(0.5, 0.0025, txt1, 'Color', 'blue');
text(0.4, 0.004, txt2, 'Color', 'red');
text(0.8, 0.008, txt3, 'Color', 'black');
text(0.8, 0.0087, txt4, 'Color', 'black');
text(0.8, 0.0013, txt5, 'Color', 'green');
% 设置标签和网格
xlabel('<P>');
ylabel('\deltaA_{ LR }');
%grid on;

% 设置纵坐标和横坐标范围
ylim([0.001 0.01]);
xlim([0 1]);

hold on;

