clear;

% 定义 alr 的函数
f_alr = @(s1,s2,s3,s4)sqrt(((s1+s2-s3-s4).*(-s1+s2-s3+s4))./((s1+s2+s3+s4).*(-s1+s2+s3-s4)));
% 初始化参数
P_plus_values = linspace(0.1, 1, 100);  % P_plus 的值
p_minu_values = linspace(0.1, 1, 100);  % p_minu 的值
ALR = 0.15;
total_N = 10^9;

% 初始化结果存储
delta_A_LR_values_optimized = zeros(length(P_plus_values), length(p_minu_values));
delta_A_LR_value_optimized_ga = zeros(length(P_plus_values), length(p_minu_values));
N1_opt_values = zeros(length(P_plus_values), length(p_minu_values));
N2_opt_values = zeros(length(P_plus_values), length(p_minu_values));
N3_opt_values = zeros(length(P_plus_values), length(p_minu_values));
N4_opt_values = zeros(length(P_plus_values), length(p_minu_values));

% 开始迭代计算
for i = 1:length(P_plus_values)
    P_plus = P_plus_values(i);
    for j = 1:length(p_minu_values)
        p_minu = p_minu_values(j);

        % 定义目标函数
        objective = @(N) sqrt( ...
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

        % 使用 GA 优化
        [N_opt_ga, delta_A_LR_optimized_ga] = ga(objective, 4, [], [], Aeq, beq, lb, [], [], options_ga);

        % 储存 GA 优化结果
        delta_A_LR_value_optimized_ga(i, j) = delta_A_LR_optimized_ga;

        % 设置 fmincon 的优化选项
        options_fmincon = optimoptions('fmincon', ...
            'Algorithm', 'sqp', ...  % 使用顺序二次规划（SQP）算法
            'Display', 'iter', ...  % 显示迭代信息
            'UseParallel', true);  % 使用并行计算

        % 使用 fmincon 进一步优化
        [N_opt, delta_A_LR_optimized] = fmincon(objective, N_opt_ga, [], [], Aeq, beq, lb, [], [], options_fmincon);

        % 存储 fmincon 优化的计算结果
        delta_A_LR_values_optimized(i, j) = delta_A_LR_optimized;
        N1_opt_values(i, j) = N_opt(1);
        N2_opt_values(i, j) = N_opt(2);
        N3_opt_values(i, j) = N_opt(3);
        N4_opt_values(i, j) = N_opt(4);
    end
end

% 结果绘图部分保持不变
[X, Y] = meshgrid(P_plus_values, p_minu_values);

% 绘制 delta_A_LR 优化值的三维图
figure;
surf(X, Y, delta_A_LR_values_optimized);
xlabel('$P^{+}$', 'Interpreter', 'latex');
ylabel('$P^{-}$', 'Interpreter', 'latex');
zlabel('$\delta A_{LR}$', 'Interpreter', 'latex');
title('Optimized \delta A_{LR}', 'Interpreter', 'latex');

% 绘制优化分配的 N1, N2, N3, N4
figure;
subplot(2, 2, 1);
surf(X, Y, N1_opt_values);
title('Optimized N1', 'Interpreter', 'latex');
xlabel('$P^{+}$', 'Interpreter', 'latex');
ylabel('$P^{-}$', 'Interpreter', 'latex');
zlabel('N1', 'Interpreter', 'latex');

subplot(2, 2, 2);
surf(X, Y, N2_opt_values);
title('Optimized N2', 'Interpreter', 'latex');
xlabel('$P^{+}$', 'Interpreter', 'latex');
ylabel('$P^{-}$', 'Interpreter', 'latex');
zlabel('N2', 'Interpreter', 'latex');

subplot(2, 2, 3);
surf(X, Y, N3_opt_values);
title('Optimized N3', 'Interpreter', 'latex');
xlabel('$P^{+}$', 'Interpreter', 'latex');
ylabel('$P^{-}$', 'Interpreter', 'latex');
zlabel('N3', 'Interpreter', 'latex');

subplot(2, 2, 4);
surf(X, Y, N4_opt_values);
title('Optimized N4', 'Interpreter', 'latex');
xlabel('$P^{+}$', 'Interpreter', 'latex');
ylabel('$P^{-}$', 'Interpreter', 'latex');
zlabel('N4', 'Interpreter', 'latex');
