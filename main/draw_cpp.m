clear all;
close all;
data1=load('./data/33_100_control.txt');
data2=load('./data/33_100_vehicle.txt'); 
% data1=load('result_control2.txt');
% data2=load('result_vehicle2.txt'); 

num = 2;

col = size(data1, 1);

Delta_1=data1(:,1);     Delta_2=data1(:,2);     Delta_3=data1(:,3);
s_Alpha = data1(:,4);   s_Beta = data1(:,5);    s_Mu = data1(:,6);
M_x = data1(:, 7);      M_y = data1(:, 8);      M_z = data1(:, 9);

T=data2(:,1);           V=data2(:,2);           Ma=data2(:,3);
Alpha=data2(:,4);       Beta=data2(:,5);        Mu=data2(:,6);
X=data2(:,7);           Y=data2(:,8);           Z=data2(:,9);
Alpha_ref=data2(:,10);   Beta_ref=data2(:,11);    Mu_ref=data2(:,12);

% index = find(abs(Alpha - 3) >= 0.05 * abs(3 - Alpha(1)), 1, 'last');
% settling_time = T(index);
% disp(['调节时间为：', num2str(settling_time), ' 秒']);

figure;
subplot(3, 1, 1);
plot(T, Alpha, 'r', T, Alpha_ref, 'b--', 'LineWidth', 1.5);
ylabel('\alpha(deg)');
legend('输出姿态角', '期望姿态角');
subplot(3, 1, 2);
plot(T, Beta, 'r', T, Beta_ref, 'b--',  'LineWidth', 1.5);
ylabel('\beta(deg)');
subplot(3, 1, 3);
plot(T, Mu, 'r', T, Mu_ref, 'b--',  'LineWidth', 1.5);
ylabel('\mu(deg)');
xlabel('t/s');

if num == 2
    axes('position',[0.5 0.75 0.075 0.075]);
    plot(T, Alpha, 'r', T, Alpha_ref, 'b--', 'LineWidth', 1.5);
    xlim([9.5,12]);

    axes('position',[0.5 0.52 0.075 0.075]);
    plot(T, Beta, 'r', T, Beta_ref, 'b--',  'LineWidth', 1.5);
    xlim([9.5,12]);
end

% 
if num == 2
    figure;
    subplot(3, 1, 1);
    plot(T, s_Alpha, 'r', 'LineWidth', 1.5);
    ylabel('s(\alpha)');
    legend('姿态角滑模函数');
    subplot(3, 1, 2);
    plot(T, s_Beta, 'r',  'LineWidth', 1.5);
    ylabel('s(\beta)');
    subplot(3, 1, 3);
    plot(T, s_Mu, 'r',  'LineWidth', 1.5);
    ylabel('s(\mu)');
    xlabel('t/s');
    axis
end

figure;
subplot(3, 1, 1);
plot(T, Delta_1, 'r', 'LineWidth', 1.5);
ylabel('\delta_e(deg)');
legend('舵面偏角');
subplot(3, 1, 2);
plot(T, Delta_2, 'r',  'LineWidth', 1.5);
ylabel('\delta_a(deg)');
subplot(3, 1, 3);
plot(T, Delta_3, 'r',  'LineWidth', 1.5);
ylabel('\delta_r(deg)');
xlabel('t/s');
% 
if num == 2
    figure;
    subplot(3, 1, 1);
    plot(T, M_x, 'r', 'LineWidth', 1.5);
    ylabel('M_x(N·m)');
    legend('实际力矩');
    subplot(3, 1, 2);
    plot(T, M_y, 'r',  'LineWidth', 1.5);
    ylabel('M_y(N·m)');
    subplot(3, 1, 3);
    plot(T, M_z, 'r',  'LineWidth', 1.5);
    ylabel('M_z(N·m)');
    xlabel('t/s');
end


% figure;
% plot(x_T,y_X,'r--');
% legend('X');
% figure;
% plot(x_T,y_Y,'r--');
% legend('Y');
% figure;
% plot(x_T,y_Z,'r--');
% legend('Z');
% figure;
% plot(x_T,y_V,'r--');
% legend('V');
% 
% data3 = load('result_control1.txt');
% data4 = load('result_vehicle1.txt'); 
% Alpha1 = data4(:,4);
% Beta1 = data4(:,5);
% Mu1 = data4(:,6);
% 
% data5 = load('result_control2.txt');
% data6 = load('result_vehicle2.txt'); 
% Alpha2 = data6(:,4);
% Beta2 = data6(:,5);
% Mu2 = data6(:,6);

% figure();
% subplot(3, 1, 1);
% plot(T, Alpha, 'r', T, Alpha1, 'b', T, Alpha2, 'm', T, Alpha_Ref, 'k--',  'LineWidth', 1.5);
% ylabel('α(deg)');
% pos = legend('基于FPSO的AC算法参数调节姿态角','参数固定姿态角','AC算法参数调节姿态角');
% set(pos,'position',[0.55,0.73,0.3,0.1]);
% subplot(3, 1, 2);
% plot(T, Beta, 'r', T, Beta1, 'b', T, Beta2, 'm', T, Beta_Ref, 'k--', 'LineWidth', 1.5);
% ylabel('β(deg)');
% subplot(3, 1, 3);
% plot(T, Mu, 'r', T, Mu1, 'b',T, Mu2, 'm', T, Mu_Ref, 'k--', 'LineWidth', 1.5);
% ylabel('γ_c(deg)');
% xlabel('t/s');
% 
% axes('position',[0.4 0.75 0.075 0.075]);
% plot(T, Alpha, 'r', T, Alpha1, 'b', T, Alpha2, 'm', T, Alpha_Ref, 'k--',  'LineWidth', 1.5);
% xlim([9.5,12]);
% 
% axes('position',[0.5 0.52 0.075 0.075]);
% plot(T, Beta, 'r', T, Beta1, 'b', T, Beta2, 'm', T, Beta_Ref, 'k--', 'LineWidth', 1.5);
% xlim([9.5,12]);

