clear;
syms dp P_minu N dp A_LR N3 N4
% 变量赋值
A_LR_val = 0.16;
N_val = 1 * 3 * 10^6;
N_val4=0.5*3*10^6;
N3_val4=0.25*3*10^6;
N4_val4=0.25*3*10^6;
P_minu_val = linspace(0.1, 0.8, 1000);

%定义函数
delta_alr=sqrt(1/((P_minu^2)*N)+(A_LR*dp)^2);
delta_alrf=sqrt(1/((P_minu^2)*N)+(A_LR*(((1-P_minu^2)/(2*P_minu^2))*sqrt(1/N3+1/N4)))^2);

num_delta_alrf=matlabFunction(delta_alrf,'Var',{P_minu,A_LR,N,N3,N4});
num_delta_alr=matlabFunction(delta_alr,'Var',{P_minu,dp,A_LR,N});


%计算具体值
num_delta_alr4=num_delta_alrf(P_minu_val,A_LR_val,N_val4,N3_val4,N4_val4);

dp_val1=0.01;
num_delta_alr1=num_delta_alr(P_minu_val,dp_val1,A_LR_val,N_val);

dp_val1=0.02;
num_delta_alr2=num_delta_alr(P_minu_val,dp_val1,A_LR_val,N_val);

figure;
plot(P_minu_val, num_delta_alr4, 'g-');
hold on;
plot(P_minu_val, num_delta_alr1, 'b-', 'DisplayName', 'dp = 0.01');
hold on;
plot(P_minu_val, num_delta_alr2, 'r--', 'DisplayName', 'dp = 0.02');
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

% 设置纵坐标和横坐标范围
ylim([0.001 0.01]);
xlim([0 1]);

