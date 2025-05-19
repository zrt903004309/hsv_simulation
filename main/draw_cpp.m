clear all;
close all;
data1=load('./data/35_100_control.txt');
data2=load('./data/35_100_vehicle.txt'); 
% data1=load('result_control2.txt');
% data2=load('result_vehicle2.txt'); 

num = 2;

col = size(data1, 1);

Delta_e=data1(:,1);     Delta_a=data1(:,2);     Delta_r=data1(:,3);
s_Alpha = data1(:,4);   s_Beta = data1(:,5);    s_Mu = data1(:,6);
M_x = data1(:, 7);      M_y = data1(:, 8);      M_z = data1(:, 9);
Mc_xe = data1(:, 10);   Mc_ye = data1(:, 11);   Mc_ze = data1(:, 12);
Mc_x = data1(:, 13);    Mc_y = data1(:, 14);    Mc_z = data1(:, 15);
u1 = data1(:, 16);      u2 = data1(:, 17);      u3 = data1(:, 18); 
% A_K = data1(:, 19);     B_K = data1(:, 20);     M_K = data1(:, 21);

T=data2(:,1);           V=data2(:,2);           Ma=data2(:,3);
Alpha=data2(:,4);       Beta=data2(:,5);        Mu=data2(:,6);
X=data2(:,7);           Y=data2(:,8);           Z=data2(:,9);
Alpha_ref=data2(:,10);   Beta_ref=data2(:,11);    Mu_ref=data2(:,12);
Gamma = data2(:, 13);   Chi = data2(:, 14);
p = data2(:, 15);       q = data2(:, 16);       r = data2(:, 17);

% index = find(abs(Alpha - 3) >= 0.05 * abs(3 - Alpha(1)), 1, 'last');
% settling_time = T(index);
% disp(['调节时间为：', num2str(settling_time), ' 秒']);

% figure;
% plot(T, A_K, 'r',T, B_K, 'b',T, M_K, 'm', 'LineWidth', 1.5);
% xlabel('t/s');
% ylabel('G_1');
% legend('g_{11}', 'g_{12}', 'g_{13}');
% figure;
% subplot(3, 1, 1);
% plot(T, A_K, 'r', 'LineWidth', 1.5);
% ylabel('g_{11}');
% % legend('G1');
% subplot(3, 1, 2);
% plot(T, B_K, 'r',  'LineWidth', 1.5);
% ylabel('g_{12}');
% subplot(3, 1, 3);
% plot(T, M_K, 'r',  'LineWidth', 1.5);
% ylabel('g_{13}');
% xlabel('t/s');
% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter4\G1变化曲线');

% fig1 - 三通道姿态角与参考姿态角
figure;
% set(gca,'FontSize',20); %设置坐标轴字体大小为8
% set(legend,'FontSize',8);   %设置legend的字体大小为8
subplot(3, 1, 1);
plot(T, Alpha, 'r', T, Alpha_ref, 'b--', 'LineWidth', 1.5);
ylim([-3,12]);
ylabel('\alpha(°)');
legend('输出姿态角', '期望姿态角','Location','southeast');
subplot(3, 1, 2);
plot(T, Beta, 'r', T, Beta_ref, 'b--',  'LineWidth', 1.5);
ylabel('\beta(°)');
legend('输出姿态角', '期望姿态角');
subplot(3, 1, 3);
plot(T, Mu, 'r', T, Mu_ref, 'b--',  'LineWidth', 1.5);
ylabel('\mu(°)');
legend('输出姿态角', '期望姿态角');
xlabel('t/s');

if num == 2
    axes('position',[0.5 0.75 0.075 0.075]);
    plot(T, Alpha, 'r', T, Alpha_ref, 'b--', 'LineWidth', 1.5);
    xlim([14,17]);

    axes('position',[0.5 0.52 0.075 0.075]);
    plot(T, Beta, 'r', T, Beta_ref, 'b--',  'LineWidth', 1.5);
    xlim([14,17]);
    
    axes('position',[0.5 0.22 0.075 0.075]);
    plot(T, Mu, 'r', T, Mu_ref, 'b--',  'LineWidth', 1.5);
    xlim([14,17]);
end

% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter6\故障姿态角跟踪响应曲线');

% fig2 - 三通道滑模函数
figure;
subplot(3, 1, 1);
plot(T, s_Alpha, 'r', 'LineWidth', 1.5);
ylim([-0.05,0.05]);
ylabel('s(\alpha)');
% legend('姿态角滑模函数');
subplot(3, 1, 2);
plot(T, s_Beta, 'r',  'LineWidth', 1.5);
ylim([-0.05,0.05]);
ylabel('s(\beta)');
subplot(3, 1, 3);
plot(T, s_Mu, 'r',  'LineWidth', 1.5);
ylim([-0.05,0.05]);
ylabel('s(\mu)');
xlabel('t/s');
axis
% axes('position',[0.45 0.84 0.15 0.09]);
% plot(T, s_Alpha, 'r', 'LineWidth', 1.5);
% xlim([14,17]);
% 
% axes('position',[0.45 0.53 0.15 0.09]);
% plot(T, s_Beta, 'r',  'LineWidth', 1.5);
% xlim([14,17]);
% 
% axes('position',[0.45 0.15 0.15 0.1]);
% plot(T, s_Mu, 'r',  'LineWidth', 1.5);
% xlim([14,17]);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter6\故障滑模函数曲线');

% fig3 - 三个舵偏角
figure;
subplot(3, 1, 1);
plot(T, Delta_e, 'r', 'LineWidth', 1.5);
ylabel('\delta_e(°)');
% legend('舵面偏角');
subplot(3, 1, 2);
plot(T, Delta_a, 'r',  'LineWidth', 1.5);
ylabel('\delta_a(°)');
subplot(3, 1, 3);
plot(T, Delta_r, 'r',  'LineWidth', 1.5);
ylabel('\delta_r(°)');
xlabel('t/s');

axes('position',[0.45 0.84 0.15 0.09]);
plot(T, Delta_e, 'r', 'LineWidth', 1.5);
xlim([14,17]);

axes('position',[0.45 0.53 0.15 0.09]);
plot(T, Delta_a, 'r',  'LineWidth', 1.5);
xlim([14,17]);

axes('position',[0.45 0.15 0.15 0.1]);
plot(T, Delta_r, 'r',  'LineWidth', 1.5);
xlim([14,17]);

% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter6\故障舵偏角曲线');

% fig4 -期望实际力矩曲线
figure;
subplot(3, 1, 1);
plot(T, M_x, 'r', 'LineWidth', 1.5);
ylabel('M_x(N·m)');
legend('期望实际力矩');
subplot(3, 1, 2);
plot(T, M_y, 'r',  'LineWidth', 1.5);
ylabel('M_y(N·m)');
subplot(3, 1, 3);
plot(T, M_z, 'r',  'LineWidth', 1.5);
ylabel('M_z(N·m)');
xlabel('t/s');

% fig5 - 控制力矩曲线
figure;
subplot(3, 1, 1);
plot(T, Mc_xe, 'r', T, Mc_x, 'b--', 'LineWidth', 1.5);
ylabel('Mc_x(N·m)');
legend('期望控制力矩', '实际控制力矩');
subplot(3, 1, 2);
plot(T, Mc_ye, 'm', T, Mc_y, 'b--', 'LineWidth', 1.5);
legend('期望控制力矩', '实际控制力矩');
ylabel('Mc_y(N·m)');
subplot(3, 1, 3);
plot(T, Mc_ze, 'r', T, Mc_z, 'b--',  'LineWidth', 1.5);
legend('期望控制力矩', '实际控制力矩');
ylabel('Mc_z(N·m)');
xlabel('t/s');


% fig6 - 速度变化曲线
figure;
plot(T,V,'r--');
legend('V');

% fig7 - 三轴位置变化曲线
figure;
plot3(X, Y, -Z, 'b-', 'LineWidth', 0.8, 'MarkerSize', 8);

% 设置坐标轴标签和标题
xlabel('X');
ylabel('Y');
zlabel('Z');
title('三维曲线图');

% 显示网格线
grid on;

% fig8 - p q r曲线
figure;
subplot(3, 1, 1);
plot(T, p, 'r', 'LineWidth', 1.5);
ylabel('p(rad/s)');
% legend('角加速度');
subplot(3, 1, 2);
plot(T, q, 'r',  'LineWidth', 1.5);
ylabel('q(rad/s)');
subplot(3, 1, 3);
plot(T, r, 'r',  'LineWidth', 1.5);
ylabel('r(rad/s)');
xlabel('t/s');
% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter6\故障姿态角速度变化曲线');

% fig9 - 航迹角
figure;
plot(T,Gamma,'r-',T,Chi,'m-');
legend('Gamma','Chi');grid on;

% fig10 - 三通道姿态角误差曲线
figure;
plot(T, Alpha-Alpha_ref, 'r',T, Beta-Beta_ref, 'b',T, Mu-Mu_ref, 'm', 'LineWidth', 1.5);
xlabel('t/s');
ylabel('姿态角误差(°)');
legend('e_\alpha', 'e_\beta', 'e_\mu');
axes('position',[0.4 0.65 0.25 0.25]);
plot(T, Alpha-Alpha_ref, 'r', T, Beta-Beta_ref, 'b--', T, Mu-Mu_ref, 'm', 'LineWidth', 1.5);
xlim([12.5,17.5]);

% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter6\故障姿态角误差');

% fig11 - 虚拟控制量曲线
% figure;
% subplot(3, 1, 1);
% plot(T, u1, 'r', 'LineWidth', 1.5);
% ylabel('u1');
% legend('期望实际力矩');
% subplot(3, 1, 2);
% plot(T, u2, 'r',  'LineWidth', 1.5);
% ylabel('u2');
% subplot(3, 1, 3);
% plot(T, u3, 'r',  'LineWidth', 1.5);
% ylabel('u3');
% xlabel('t/s');

mfd1 = data1(:, 19);    mfd2 = data1(:, 20);    mfd3 = data1(:, 21);
pfd1 = data1(:, 22);    pfd2 = data1(:, 23);    pfd3 = data1(:, 24);
figure;
plot(T, mfd1, 'r',T, mfd2, 'b',T, mfd3, 'm', 'LineWidth', 1.5);
xlabel('t/s');
ylabel('测量干扰');
ylim([-100,100]);
legend('md1', 'md2', 'md3');

figure;
plot(T, pfd1/30, 'r',T, pfd2/30, 'b',T, pfd3/30, 'm', 'LineWidth', 1.5);
xlabel('t/s');
ylabel('$\hat{f_{d}}$', 'Interpreter','latex');
ylim([-1,1]);
legend('$\hat{f_{d1}}$', '$\hat{f_{d2}}$', '$\hat{f_{d3}}$', 'Interpreter','latex');
set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter5\35km预测干扰曲线');

figure;
plot(T, pfd1-mfd1, 'r',T, pfd2-mfd2, 'b',T, pfd3-mfd3, 'm', 'LineWidth', 1.5);
xlabel('t/s');
ylabel('预测-测量干扰');
ylim([-100,100]);
legend('dd1', 'dd2', 'dd3');

% A_eta = data1(:, 25);    B_eta = data1(:, 26);    M_eta = data1(:, 27);
% A_eta = data1(:, 22);    B_eta = data1(:, 23);    M_eta = data1(:, 24);
% figure;
% plot(T, A_eta, 'r',T, B_eta, 'b',T, M_eta, 'm', 'LineWidth', 1.5);
% xlabel('t/s');
% ylabel('$\hat{\eta}$', 'Interpreter','latex');
% % ylim([-100,100]);
% legend('$\hat{\eta_1}$', '$\hat{\eta_2}$', '$\hat{\eta_3}$', 'Interpreter','latex');
% set(gcf,'PaperUnits','centimeters','PaperPosition',[14 19 14 8.6],'PaperPositionMode', 'manual');
% print(gcf,'-dpng','-r300','F:\MasterEssay\essay\thesis\figures\chapter6\切换增益变化曲线');

% data3 = load('./data/useful/三四章/pid+dann35_100_control.txt');
% data4 = load('./data/useful/三四章/pid+dann35_100_vehicle.txt'); 
% Alpha1 = data4(:,4);
% Beta1 = data4(:,5);
% Mu1 = data4(:,6);
% 
% data5 = load('./data/useful/三四章/npid+dann35_100_control.txt');
% data6 = load('./data/useful/三四章/npid+dann35_100_vehicle.txt'); 
% Alpha2 = data6(:,4);
% Beta2 = data6(:,5);
% Mu2 = data6(:,6);
% 
% figure();
% subplot(3, 1, 1);
% plot(T, Alpha, 'r', T, Alpha1, 'b', T, Alpha2, 'm', T, Alpha_ref, 'k--',  'LineWidth', 1.5);
% ylabel('α(°)');
% pos = legend('基于FPSO的AC算法参数调节姿态角','参数固定姿态角','AC算法参数调节姿态角');
% set(pos,'position',[0.55,0.73,0.3,0.1]);
% subplot(3, 1, 2);
% plot(T, Beta, 'r', T, Beta1, 'b', T, Beta2, 'm', T, Beta_ref, 'k--', 'LineWidth', 1.5);
% ylabel('β(°)');
% subplot(3, 1, 3);
% plot(T, Mu, 'r', T, Mu1, 'b',T, Mu2, 'm', T, Mu_ref, 'k--', 'LineWidth', 1.5);
% ylabel('γ_c(°)');
% xlabel('t/s');
% 
% axes('position',[0.4 0.75 0.075 0.075]);
% plot(T, Alpha, 'r', T, Alpha1, 'b', T, Alpha2, 'm', T, Alpha_Ref, 'k--',  'LineWidth', 1.5);
% xlim([9.5,12]);
% 
% axes('position',[0.5 0.52 0.075 0.075]);
% plot(T, Beta, 'r', T, Beta1, 'b', T, Beta2, 'm', T, Beta_Ref, 'k--', 'LineWidth', 1.5);
% xlim([9.5,12]);

