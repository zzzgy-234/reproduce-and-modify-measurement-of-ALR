% 定义ALR关于四个散射截面的函数用于误差传递的方法
f_alr = @(s1, s2, s3, s4) (4 .* (s1 - s3) .* (s2 - s3)) ./ (s3 .* (s1 - s2)) - ...
    (4 .* (s1 - s3) .* (s2 - s3) .* ((s1 .* (-s3 ./ (4 .* (s1 - s3) .* (s2 - s3) .* (s1 + s2 - s3 - s4))).^(1/2)) ./ 2 - ...
    (s2 .* (-s3 ./ (4 .* (s1 - s3) .* (s2 - s3) .* (s1 + s2 - s3 - s4))).^(1/2)) ./ 2 + 1)) ./ (s3 .* (s1 - s2));

% 定义ALR关于P和散射截面的函数
f_alr1 = @(P_minu, sa, sb) (1./P_minu) .* ((sb - sa) ./ (sa + sb));

% 初始化数组长度
P = linspace(0.1, 0.8, 1000); % 生成从0.1到0.8的10个点
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
    lambda_na = 0.5 * 3 * (10^6); % na的平均值
    lambda_nb = 0.5 * 3 * (10^6); % nb的平均值

    num_samples = 100000; % 生成的样本数量

    % 生成泊松分布的随机样本
    n1_sample = poissrnd(lambda_n1, num_samples, 1);
    n2_sample = poissrnd(lambda_n2, num_samples, 1);
    n3_sample = poissrnd(lambda_n3, num_samples, 1);
    n4_sample = poissrnd(lambda_n4, num_samples, 1);
    na_sample = poissrnd(lambda_na, num_samples, 1);
    nb_sample = poissrnd(lambda_nb, num_samples, 1);

    % 生成高斯分布的随机样本
    sigma1 = 0.01 * P(i); % 标准差
    sigma2 = 0.02 * P(i); % 标准差
    P1_minu_samples = P(i) + sigma1 * randn(num_samples, 1);
    P2_minu_samples = P(i) + sigma2 * randn(num_samples, 1);

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

    % 计算平均值和标准差
    f_alr_mean(i) = mean(f1_samples);
    f_alr_std_erro(i) = std(f1_samples);
    f_alr1_mean(i) = mean(f2_samples);
    f_alr1_std_erro(i) = std(f2_samples);
    f_alr2_mean(i) = mean(f3_samples);
    f_alr2_std_erro(i) = std(f3_samples);
end

% 绘制误差曲线
figure;
plot(P, f_alr_std_erro, 'r', 'DisplayName', 'Error for f_{ALR}'); % 红色
xlabel('P');
ylabel('\deltaA_{LR}');
hold on;
plot(P, f_alr1_std_erro, 'g', 'DisplayName', 'Error for f_{ALR1}'); % 绿色
plot(P, f_alr2_std_erro, 'b', 'DisplayName', 'Error for f_{ALR2}'); % 蓝色
xlim([0 1]); % 设置横坐标范围
ylim([0 0.01]); % 设置纵坐标范围
legend show; % 显示图例
grid on;
hold off;

% 绘制A_{LR}真值以及模拟值
figure;
% 绘制真值 0.16 作为参考线
plot(P, 0.16 * ones(size(P)), 'k--', 'DisplayName', 'True Value A_{LR} = 0.16'); % 黑色虚线表示真值
xlabel('P');
ylabel('A_{LR}');
hold on;

% 绘制误差包络线：真值 0.16 加上/减去误差
plot(P, f_alr_mean, 'r', 'DisplayName', 'f_{ALR1}'); % 红色
fill([P, fliplr(P)], [0.16 + f_alr_std_erro, fliplr(0.16 - f_alr_std_erro)], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 真值0.16的误差包络线

plot(P, f_alr1_mean, 'g', 'DisplayName', 'f_{ALR1}'); % 绿色
fill([P, fliplr(P)], [0.16 + f_alr1_std_erro, fliplr(0.16 - f_alr1_std_erro)], ...
    'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % f_{ALR1}的误差包络线

plot(P, f_alr2_mean, 'b', 'DisplayName', 'f_{ALR2}'); % 蓝色
fill([P, fliplr(P)], [0.16 + f_alr2_std_erro, fliplr(0.16 - f_alr2_std_erro)], ...
    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % f_{ALR2}的误差包络线

legend show; % 显示图例
grid on;
hold off;
