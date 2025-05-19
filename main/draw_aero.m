clear all;
close all;
data=load('./data/aero_coefficients.txt');

Ma = data(:, 1);    Alpha = data(:, 2);
CD = data(:, 3);    CL = data(:, 4);    Cm_alpha = data(:, 5);  LD_ratio = data(:, 6);
Cl_LE = data(:, 7); Cm_LE = data(:, 8); Cn_LE = data(:, 9);
Cl_RUD = data(:, 10);Cm_RUD = data(:, 11);Cn_RUD = data(:, 12);

% 确定x和y的唯一值
unique_Ma = unique(Ma);
unique_Alpha = unique(Alpha);

% 创建网格
[ma, alpha] = meshgrid(unique_Ma, unique_Alpha);

% 初始化系数矩阵
cd = nan(size(ma)); cl = nan(size(ma)); cm_a = nan(size(ma)); ld_ratio = nan(size(ma));
cl_le = nan(size(ma));cm_le = nan(size(ma));cn_le = nan(size(ma));
cl_rud = nan(size(ma));cm_rud = nan(size(ma));cn_rud = nan(size(ma));

% 填充Z矩阵
for i = 1:length(unique_Ma)
    for j = 1:length(unique_Alpha)
        index = (i-1)*length(unique_Alpha) + j;
        cd(j, i) = CD(index);   cl(j, i) = CL(index);   cm_a(j, i) = Cm_alpha(index);   ld_ratio(j, i) = LD_ratio(index);
        cl_le(j, i) = Cl_LE(index); cm_le(j, i) = Cm_LE(index); cn_le(j, i) = Cn_LE(index);
        cl_rud(j, i) = Cl_RUD(index);   cm_rud(j, i) = Cm_RUD(index);   cn_rud(j, i) = Cn_RUD(index);
    end
end

Fig1 = mesh(alpha, ma, cd);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('CD_\alpha', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig1,'./fig/chapter2/CD_alpha.png');

Fig2 = mesh(alpha, ma, cl);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('CL_\alpha', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig2,'./fig/chapter2/CL_alpha.png');

Fig3 = mesh(alpha, ma, cm_a);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cm_\alpha', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig3,'./fig/chapter2/Cm_alpha.png');

Fig4 = mesh(alpha, ma, cl_le);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cl-\delta_e', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig4,'./fig/chapter2/Cl_deltae.png');

Fig5 = mesh(alpha, ma, cm_le);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cm-\delta_e', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig5,'./fig/chapter2/Cm_deltae.png');

Fig6 = mesh(alpha, ma, cn_le);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cn-\delta_e', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig6,'./fig/chapter2/Cn_deltae.png');

Fig7 = mesh(alpha, ma, cl_rud);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cl-\delta_r', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig7,'./fig/chapter2/Cl_deltar.png');

Fig8 = mesh(alpha, ma, cm_rud);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cm-\delta_r', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig8,'./fig/chapter2/Cm_deltar.png');

Fig9 = mesh(alpha, ma, cn_rud);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('Cn-\delta_r', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig9,'./fig/chapter2/Cn_deltar.png');

Fig10 = mesh(alpha, ma, ld_ratio);
xlabel('\alpha(deg)', 'FontSize',18);
ylabel('Ma', 'FontSize',18);
zlabel('L-D ratio', 'FontSize',18);

set(gca, 'FontSize', 14);
% set(gca, 'LineWidth', 0.8); % 设置刻度线长度
view(-40, 10);

saveas(Fig10,'./fig/chapter2/LD_ratio.png');