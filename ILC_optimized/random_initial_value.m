clear;
%定义alr的函数
f_alr = @(s1,s2,s3,s4)sqrt(((s1+s2-s3-s4).*(-s1+s2-s3+s4))./((s1+s2+s3+s4).*(-s1+s2+s3-s4)));
p_minu = 0.8;
ALR = 0.15;
P_plus_values = linspace(0.1, 1, 10);
alr_mean = zeros(1,10);
alr_std_erro = zeros(1,10);
total_N = 10^10;
% 在循环外定义一个数组来存储每次循环的 N1_opt_values
all_N1_values = zeros(5, length(P_plus_values));
all_delta_alr_values = zeros(5, length(P_plus_values));
for repeat = 1:5
    fprintf('Iteration: %d\n', repeat);
    
    alr_mean = zeros(1,10);
    alr_std_erro = zeros(1,10);
    
    % 初始化存储 delta_A_LR 和 N 值的数组
    delta_A_LR_values_optimized = zeros(size(P_plus_values));
    delta_A_LR_values_unoptimized = zeros(size(P_plus_values));
    delta_A_LR_value_optimized_ga = zeros(size(P_plus_values));
    N1_opt_values = zeros(size(P_plus_values));
    N2_opt_values = zeros(size(P_plus_values));
    N3_opt_values = zeros(size(P_plus_values));
    N4_opt_values = zeros(size(P_plus_values));

    for i = 1:length(P_plus_values)
        P_plus = P_plus_values(i);

    % 随机扰动均匀分配的 N1, N2, N3, N4
    random_disturbance = rand(1, 4);  % 生成 4 个随机数
    random_disturbance = random_disturbance / sum(random_disturbance);  % 归一化，使其和为 1
    N_unoptimized = total_N * random_disturbance;  % 将扰动应用到总量上4];

        % 计算未优化的 delta_A_LR (均匀分配 N 的情况)
    delta_1 = (1 ./ (8 .* P_plus .* p_minu)) .* (-P_plus + p_minu - (P_plus .* p_minu - 1) .* ALR) .* (1 - P_plus .* p_minu + ALR .* (P_plus - p_minu)) .* (1 ./ sqrt(N_unoptimized(1)));
    delta_2 = (1 ./ (8 .* P_plus .* p_minu)) .* (-P_plus - p_minu - (P_plus .* p_minu + 1) .* ALR) .* (1 + P_plus .* p_minu + ALR .* (-P_plus - p_minu)) .* (1 ./ sqrt(N_unoptimized(2)));
    delta_3 = (1 ./ (8 .* P_plus .* p_minu)) .* (P_plus + p_minu - (P_plus .* p_minu + 1) .* ALR) .* (1 + P_plus .* p_minu + ALR .* (P_plus + p_minu)) .* (1 ./ sqrt(N_unoptimized(3)));
    delta_4 = (1 ./ (8 .* P_plus .* p_minu)) .* (P_plus - p_minu - (P_plus .* p_minu - 1) .* ALR) .* (1 - P_plus .* p_minu + ALR .* (-P_plus + p_minu)) .* (1 ./ sqrt(N_unoptimized(4)));
    delta_A_LR_unoptimized = sqrt(delta_1.^2 + delta_2.^2 + delta_3.^2 + delta_4.^2) * 1e5;
    
    % 存储未优化的计算结果
    delta_A_LR_values_unoptimized(i) = delta_A_LR_unoptimized;
    
    % 定义目标函数
    objective = @(N) sqrt(...
        ((1 ./ (8 .* P_plus .* p_minu)) .* (-P_plus + p_minu - (P_plus .* p_minu - 1) .* ALR) .* (1 - P_plus .* p_minu + ALR .* (P_plus - p_minu)) .* (1 ./ sqrt(N(1))))^2 + ...
        ((1 ./ (8 .* P_plus .* p_minu)) .* (-P_plus - p_minu - (P_plus .* p_minu + 1) .* ALR) .* (1 + P_plus .* p_minu + ALR .* (-P_plus - p_minu)) .* (1 ./ sqrt(N(2))))^2 + ...
        ((1 ./ (8 .* P_plus .* p_minu)) .* (P_plus + p_minu - (P_plus .* p_minu + 1) .* ALR) .* (1 + P_plus .* p_minu + ALR .* (P_plus + p_minu)) .* (1 ./ sqrt(N(3))))^2 + ...
        ((1 ./ (8 .* P_plus .* p_minu)) .* (P_plus - p_minu - (P_plus .* p_minu - 1) .* ALR) .* (1 - P_plus .* p_minu + ALR .* (-P_plus + p_minu)) .* (1 ./ sqrt(N(4))))^2 ...
    ) * 1e5;

    % 初始值
    N0 = [2*10^8, 10^8, 10^8, 10^8];
    
    % 线性约束：N1 + N2 + N3 + N4 = 10^9
    Aeq = [1, 1, 1, 1];
    beq = 10^9;
    
    % 下界：N1, N2, N3, N4 都大于一个较小的正数
    lb = [1, 1, 1, 1];
    
% 增加遗传算法的种群规模和代数
    options_ga = optimoptions('ga', ...
    'PopulationSize', 200, ...  % 增大种群规模
    'MaxGenerations', 100, ...  % 增加代数
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationadaptfeasible, 0.02}, ...  % 增加突变率
    'Display', 'iter', ...
    'UseParallel', true);  % 使用并行计算

    
    % 使用GA优化
    [N_opt_ga, delta_A_LR_optimized_ga] = ga(objective, 4, [], [], Aeq, beq, lb, [], [], options_ga);
    %储存优化结果
    delta_A_LR_value_optimized_ga(i) = delta_A_LR_optimized_ga;
    % 设置 fmincon 的优化选项
    options_fmincon = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...  % 使用顺序二次规划（SQP）算法
        'Display', 'iter', ...  % 显示迭代信息
        'UseParallel', true);  % 使用并行计算
    

        % 使用 fmincon 进行进一步优化
        [N_opt, delta_A_LR_optimized] = fmincon(objective, N_opt_ga, [], [], Aeq, beq, lb, [], [], options_fmincon);

        % 存储优化的计算结果
    delta_A_LR_values_optimized(i) = delta_A_LR_optimized;
    N1_opt_values(i) = N_opt(1);
    N2_opt_values(i) = N_opt(2);
    N3_opt_values(i) = N_opt(3);
    N4_opt_values(i) = N_opt(4);
    %定义MC模拟的事例数
    lambda_n1 = N_opt(1);
    lambda_n2 = N_opt(2);
    lambda_n3 = N_opt(3);
    lambda_n4 = N_opt(4);

    num_samples = 10^4;
    % 生成满足泊松分布的随机样本
    n1_sample = poissrnd(lambda_n1, num_samples, 1);
    n2_sample = poissrnd(lambda_n2, num_samples, 1);
    n3_sample = poissrnd(lambda_n3, num_samples, 1);
    n4_sample = poissrnd(lambda_n4, num_samples, 1);
    %计算散射截面
    s1_sample = n1_sample*(1-P_plus*p_minu+ALR*(P_plus-p_minu))/lambda_n1;
    s2_sample = n2_sample*(1+P_plus*p_minu+ALR*(-P_plus-p_minu))/lambda_n2;
    s3_sample = n3_sample*(1+P_plus*p_minu+ALR*(P_plus+p_minu))/lambda_n3;
    s4_sample = n4_sample*(1-P_plus*p_minu+ALR*(-P_plus+p_minu))/lambda_n4;
    f_alr_samples = f_alr(s1_sample,s2_sample,s3_sample,s4_sample);
    %计算
    alr_mean(i) = mean(f_alr_samples);
    alr_std_erro(i) = std(f_alr_samples)*1e5;
    % 打印调试信息
    fprintf('P_plus: %.4f, Optimized N: [%.4e, %.4e, %.4e, %.4e], delta_A_LR_optimized: %.4e\n,delta_A_LR_ga: %.4e\n, delta_A_LR_MC: %.4e\n ,alr_value : %.4e\n', ...
        P_plus, N_opt(1), N_opt(2), N_opt(3), N_opt(4), delta_A_LR_optimized,delta_A_LR_optimized_ga,alr_std_erro,alr_mean);
    end

    % 存储每次循环的 N1 值
    all_N1_values(repeat, :) = N1_opt_values;
    all_delta_alr_values(repeat, :)=delta_A_LR_values_optimized;
  
end

% 绘制所有 N1 值
figure;
hold on;
for repeat = 1:5
    plot(P_plus_values, all_N1_values(repeat, :), '-o', 'DisplayName', ['N 1Iteration ', num2str(repeat)]);
end
xlabel('P\_plus');
ylabel('N1 Values');
title('N1 Values for Each Iteration');
legend;
hold off;

% 绘制所有 delta_A_LR_value_optimized_ga 值
figure;
hold on;
for repeat = 1:5
    plot(P_plus_values, all_delta_alr_values(repeat, :), '-x', 'DisplayName', ['\Delta A_{LR} GA Iteration ', num2str(repeat)]);
end
xlabel('P\_plus');
ylabel('\Delta A_{LR} (optimized)');
title('\Delta A_{LR}Optimized Values for Each Iteration');
legend;
hold off;