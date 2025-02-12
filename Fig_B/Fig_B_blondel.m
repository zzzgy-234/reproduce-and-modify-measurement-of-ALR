clear;
clc;
syms N P_minu

% 定义方程
eqn_standard = (0.002)^2 == 1/(P_minu^2 * N) + (0.16 * 0.01)^2;
eqn_four_bunch = (0.002)^2 == 1/(P_minu^2 * (N/2)) + (0.16 * ((1 - P_minu^2) / (2 * P_minu^2)) * sqrt(1/(N/4) + 1/(N/4)))^2;

% 解方程得到N
N_standard = solve(eqn_standard, N);
N_four_bunch = solve(eqn_four_bunch, N);

% 将符号解转换为函数
N_standard_func = matlabFunction(N_standard, 'Vars', P_minu);
N_four_bunch_func = matlabFunction(N_four_bunch, 'Vars', P_minu);

% 定义P_minu的范围
P_minu_vals = linspace(0.1, 1, 100);

% 计算给定P_minu范围内的N值
N_standard_vals = N_standard_func(P_minu_vals);
N_four_bunch_vals = N_four_bunch_func(P_minu_vals);

% 绘制结果
figure;
plot(P_minu_vals, N_standard_vals, 'b', 'LineWidth', 2);
hold on;
plot(P_minu_vals, N_four_bunch_vals, 'r', 'LineWidth', 2);
set(gca, 'YScale', 'log'); % 设置纵坐标为对数尺度
xlabel('P_{minu}');
ylabel('N');
legend('standard scheme', '4 bunches scheme');

grid on;
hold off;
