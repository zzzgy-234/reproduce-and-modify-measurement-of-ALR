% 定义ALR关于四个散射截面的函数用于误差传递的方法
f_alr = @(s1, s2, s3, s4) (4 .* (s1 - s3) .* (s2 - s3)) ./ (s3 .* (s1 - s2)) - (4 .* (s1 - s3) .* (s2 - s3) .* ((s1 .* (-s3 ./ (4 .* (s1 - s3) .* (s2 - s3) .* (s1 + s2 - s3 - s4))).^(1/2)) ./ 2 - (s2 .* (-s3 ./ (4 .* (s1 - s3) .* (s2 - s3) .* (s1 + s2 - s3 - s4))).^(1/2)) ./ 2 + 1)) ./ (s3 .* (s1 - s2));
% 定义ALR关于P和散射截面的函数
f_alr1 = @(P_minu, sa, sb) (1./P_minu) .* ((sb - sa) ./ (sa + sb));

% 初始化数组长度
P = linspace(0.1, 0.8, 1000); % 生成从0.1到1的1000个点
f_alr_mean = zeros(1, 1000);
f_alr_std_erro = zeros(1, 1000);
f_alr1_mean = zeros(1, 1000);
f_alr1_std_erro = zeros(1, 1000);
f_alr2_mean = zeros(1, 1000);
f_alr2_std_erro = zeros(1, 1000);

for i = 1:1000


    % 定义泊松分布的参数
    lambda_n1 = 0.25 * 3 * (10^6); % n1的平均值
    lambda_n2 = 0.25 * 3 * (10^6); % n2的平均值
    lambda_n3 = 0.25 * 3 * (10^6); % n3的平均值
    lambda_n4 = 0.25 * 3 * (10^6); % n4的平均值
    % standard scheme 的方程泊松赋值
    lambda_na = 0.5 * 3 * (10^6);
    lambda_nb = 0.5 * 3 * (10^6);

    num_samples = 100000; 
    num_samples1 = 100000;
    % 生成的样本数量

    % 生成满足泊松分布的随机样本
    n1_sample = poissrnd(lambda_n1, num_samples, 1);
    n2_sample = poissrnd(lambda_n2, num_samples, 1);
    n3_sample = poissrnd(lambda_n3, num_samples, 1);
    n4_sample = poissrnd(lambda_n4, num_samples, 1);
    na_sample = poissrnd(lambda_na, num_samples1, 1);
    nb_sample = poissrnd(lambda_nb, num_samples1, 1);

    % 生成满足高斯分布的随机样本
    sigma1 = 0.01*P(i); % 标准差
    sigma2 = 0.02*P(i); % 标准差
    % 生成P_minu的样本点
    P1_minu_samples = P(i) + sigma1 * randn(num_samples1,1);
    P2_minu_samples = P(i) + sigma2 * randn(num_samples1,1);

    % 计算散射截面的值
    s1_samples = n1_sample .* (1 - P(i) * 0.16) / lambda_n1;
    s2_samples = n2_sample .* (1 + P(i) * 0.16) / lambda_n2;
    s3_samples = n3_sample / lambda_n3;
    s4_samples = n4_sample .* (1 - P(i) * P(i)) / lambda_n4;
    sa_samples = na_sample .* (1 - P(i) * 0.16) / lambda_na;
    sb_samples = nb_sample .* (1 + P(i) * 0.16) / lambda_nb;

    % 计算函数大小
    f1_samples = f_alr(s1_samples, s2_samples, s3_samples, s4_samples);
    f2_samples = f_alr1(P1_minu_samples, sa_samples, sb_samples);
    f3_samples = f_alr1(P2_minu_samples, sa_samples, sb_samples);

    % 计算平均值
    f_alr_mean(i) = mean(f1_samples);
    f_alr_std_erro(i) = std(f1_samples) ;
    f_alr1_mean(i) = mean(f2_samples);
    f_alr1_std_erro(i) = std(f2_samples) ;
    f_alr2_mean(i) = mean(f3_samples);
    f_alr2_std_erro(i) = std(f3_samples) ;
end


% 绘图
figure;
plot(P, f_alr_std_erro, 'r'); % 红色
xlabel('P');
ylabel('\deltaA_{LR}');
hold on;
plot(P, f_alr1_std_erro, 'g'); % 绿色
plot(P, f_alr2_std_erro, 'b'); % 蓝色
hold on;

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

