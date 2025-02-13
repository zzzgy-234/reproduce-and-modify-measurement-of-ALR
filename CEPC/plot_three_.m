

P_plus_values = linspace(0.1, 1, 100);  % P_plus 的值
p_minu_values = linspace(0.1, 1, 100);  % p_minu 的值

% 结果绘图部分保持不变
[X, Y] = meshgrid(P_plus_values, p_minu_values);

% 绘制 delta_A_LR 优化值的三维图
figure;
surf(X, Y, delta_A_LR_values_optimized*1e5);
xlabel('$P^{+}$', 'Interpreter', 'latex');
ylabel('$P^{-}$', 'Interpreter', 'latex');
zlabel('$\delta A_{LR}$', 'Interpreter', 'latex');
title('Optimized $\delta A_{LR}$', 'Interpreter', 'latex');

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