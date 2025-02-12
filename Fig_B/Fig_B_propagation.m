clear;
% 误差传递的P相对误差为1%时的值
f = @(P_minu, N, S_u, A_LR, dp) sqrt(...
    (S_u^2 * (A_LR * P_minu - 1).^2 .* (1./(P_minu .* (S_u * (A_LR * P_minu - 1) - S_u * (A_LR * P_minu + 1))) - (S_u * (A_LR * P_minu - 1) + S_u * (A_LR * P_minu + 1)) ./ (P_minu .* (S_u * (A_LR * P_minu - 1) - S_u * (A_LR * P_minu + 1))).^2).^2) ./ (0.5*N) + ...
    (S_u^2 * (A_LR * P_minu + 1).^2 .* (1./(P_minu .* (S_u * (A_LR * P_minu - 1) - S_u * (A_LR * P_minu + 1))) + (S_u * (A_LR * P_minu - 1) + S_u * (A_LR * P_minu + 1)) ./ (P_minu .* (S_u * (A_LR * P_minu - 1) - S_u * (A_LR * P_minu + 1))).^2).^2) ./(0.5*N)  + ...
    (dp^2 * (S_u * (A_LR * P_minu - 1) + S_u * (A_LR * P_minu + 1)).^2) ./ (P_minu.^2 .* (S_u * (A_LR * P_minu - 1) - S_u * (A_LR * P_minu + 1)).^2)) - 0.002;

% Define specific values for A_LR, S_u, and dp
A_LR_val = 0.16;
S_u_val = 2;
dp_val = 0.01;

% Create a function handle with fixed A_LR, S_u, and dp
f_x = @(P_minu, N) f(P_minu, N, S_u_val, A_LR_val, dp_val);

syms A_LR P_plus P_minu s1 s2 s3 s4 N1 N2 N3 N4 N % 定义参数

% 定义符号函数 f_alr
f_alr_sym = sqrt(-s1 * s2 + s1 * s3 + s2 * s3 - s3^2) / sqrt(s1 * s3 + s2 * s3 - s3^2 - s3 * s4);

% 计算 f_alr 的偏导数并平方和
grad_f_alr = gradient(f_alr_sym, [s1, s2, s3, s4]);
theo_delta_alr = sqrt((grad_f_alr(1) * sqrt(1 / N1) * s3 * (1 - P_minu * A_LR))^2 + ...
                      (grad_f_alr(2) * sqrt(1 / N2) * s3 * (1 + P_plus * A_LR))^2 + ...
                      (grad_f_alr(3) * sqrt(1 / N3) * s3)^2 + ...
                      (grad_f_alr(4) * sqrt(1 / N4) * s3 * (1 - P_plus * P_minu + (P_plus - P_minu) * A_LR))^2);

% 进行符号替换
subs_s1 = s3 * (1 - P_minu * A_LR);
subs_s2 = s3 * (1 + P_plus * A_LR);
subs_s4 = s3 * (1 - P_plus * P_minu + (P_plus - P_minu) * A_LR);
subs_N1 = 0.25 * N;
subs_N2 = 0.25* N;
subs_N3 = 0.25 * N;
subs_N4 = 0.25 * N;

theo_delta_alr = subs(theo_delta_alr, {'s1', 's2', 's4', 'N1', 'N2', 'N3', 'N4'}, {subs_s1, subs_s2, subs_s4, subs_N1, subs_N2, subs_N3, subs_N4});


% 显示新的表达式
disp('Theoretical delta for ALR:');
disp(theo_delta_alr);

% 定义参数值
% P_plus_val = 0.5;
A_LR_val = 0.16;
s3_val = 1; % 给 s3 一个具体的值
target_value = 0.002; % 目标值


theo_delta_alr=subs(theo_delta_alr, {P_plus, P_minu, s3, A_LR,N1,N2,N3,N4}, {P_plus, P_minu, s3_val, A_LR_val,1/(4*N),1/(4*N),1/(4*N),1/(4*N)});

disp('Theoretical delta for ALR:');
disp(theo_delta_alr);

eqn=[
  0.002==theo_delta_alr;
];

N_sol=solve(eqn, N);

disp(N_sol);

% 将 delta_A_LR 定义为一个函数关于 P_1
Z_events = matlabFunction(N_sol);

P_minu_val = linspace(0.1, 0.8, 1000);
P_plus_val=linspace(0.1,0.8,1000);
Z_events_val = Z_events(P_minu_val,P_plus_val);

figure;

plot(P_minu_val, Z_events_val,'b');

grid on;
hold on;
% Plot the implicit function
fimplicit(f_x, [0.1 0.8 0 10^8],'r');

xlabel('<P>');
ylabel('Zevents');
grid on;
txt1 = {'4 bunch scheme'};
txt2 = {'\deltaP/P=1%'};
txt3={'\deltaA_{LR}=0.002'};
txt4={'A_{LR}=0.16'};
text(0.4, 2*10^6, txt1, 'Color', 'blue');
text(0.44, 4*10^6, txt2, 'Color', 'red');
text(0.6, 5*10^7, txt3, 'Color', 'blue');
text(0.6, 3.6*10^7, txt4, 'Color', 'red');
% Set y-axis to logarithmic scale
set(gca, 'yscale', 'log');



