% lid driven flow算例
% 有限容积法
% 交错网格
% SIMPLE
% TODO: 1、动量求解当前还没加松弛因子，2、线性方程组求解改成迭代法，增大残差要求节约计算资源，3、当前没有检查过非稳态问题的结果
clear;

% 网格设置
n = 10;
X_domain = n;% 总长m
Y_domain = n;% 总高m
ncell_X  = n;% X向单元格数
ncell_Y  = n;% Y向单元格数

data = array_data_stagger(ncell_X, ncell_Y, X_domain, Y_domain)
% 初始化变量值
data.lam_cond = data.lam_cond*1;
data.lam_vis = data.lam_vis*0.2;

data.U_stg = (data.U_stg + 1) * 0;
data.V_stg = (data.V_stg + 1) * 0;
data.U_stg(1,:) = 1.2;
% data.P     = (data.Y);
% 时间步长
dt = 1e15; % 时间步长→∞，模拟稳态问题
i  = 1;

% load data_step1.mat

for i = 1:3e3
% 求解过程
% 更新动量方程系数
[D_U_X, D_U_Y, F_U_X, F_U_Y] = data.flux_Ugrid(data.lam_vis);%D扩导，F截面上的质量流量
[D_V_X, D_V_Y, F_V_X, F_V_Y] = data.flux_Vgrid(data.lam_vis);
[aP_U, aE_U, aW_U, aN_U, aS_U, b_U] = data.coef(D_U_X, D_U_Y, F_U_X, F_U_Y, dt, "U", data.scheme_CD);
[aP_V, aE_V, aW_V, aN_V, aS_V, b_V] = data.coef(D_V_X, D_V_Y, F_V_X, F_V_Y, dt, "V", data.scheme_CD);

% 求解动量方程
data.U_stg = data.solve_eqn(aP_U, aE_U, aW_U, aN_U, aS_U, b_U, "U");
data.V_stg = data.solve_eqn(aP_V, aE_V, aW_V, aN_V, aS_V, b_V, "V");
% 更新修正方程系数
[aP_P, aE_P, aW_P, aN_P, aS_P, b_P] = data.coef_Pprime(aP_U, aP_V);
% 锚定P'（可选，书上说会收敛更快一些）
aP_P(2,2) = 1;
aP_E(2,2) = 0;
aP_W(2,2) = 0;
aP_N(2,2) = 0;
aP_S(2,2) = 0;
b_P(2,2)  = 0;
% 求解压力修正方程
data.Pprime = data.solve_eqn(aP_P, aE_P, aW_P, aN_P, aS_P, b_P,"Pprime");
% 修正速度
[Uprime, Vprime] = data.vel_correct(aP_U, aP_V);
data.U_stg = data.U_stg + Uprime;
data.V_stg = data.V_stg + Vprime;
% 修正压力
relax_factor = 0.1;% 压力松弛因子
data.P = data.P + relax_factor .* data.Pprime;
continuity(i) = sum(abs(b_P(~isnan(b_P))));
% i = i + 1;
end

% 绘图
[UU, VV] =  data.vel_interp();
% 矢量图
figure()

colormap('summer');
contourf(data.X, data.Y, sqrt(UU.^2+VV.^2), 15,'LineStyle','none')
colorbar;
hold on
quiver(data.X, data.Y, UU, VV,"Marker",".","Color","k")
axis equal
hold off
xlim([0,n]);
ylim([0,n]);
xlabel('X/m'); % 建议加上单位
ylabel('Y/m');
% 流线图
figure()
verts = stream2(data.X, data.Y, UU, VV,data.X, data.Y,[0.1 150]);
lineobj = streamline(verts);
axis equal
xlim([0,n]);
ylim([0,n]);

% 压力云图
figure()
contourf(data.X, data.Y, data.P)
axis equal
% 压力修正值云图
figure()
contourf(data.X, data.Y, data.Pprime)
axis equal
% 连续性方程残差云图，即压力修正方程的b项
figure()
contourf(data.X, data.Y, b_P)
axis equal
% 连续性方程残差
figure()
plot(log10(continuity));
ylabel('continuity residual in log scale');

% save data_step1 data
