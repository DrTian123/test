clc
clear

mu0 = 4*pi*1e-7;  % 真空磁导率
I = 10;           % 每匝线圈电流 (A)
num_coils = 10;   % 线圈的匝数
z_target = 0;     % 中心区域高度 (m)
r_target = 0;     % 中心区域半径 (m)
target_radius = 0.01;  % 磁场聚焦区域大小 (m)

% 初始线圈参数 (匝的半径和高度)
R_init = linspace(0.01, 0.05, num_coils); % 半径初始值 (m)
Z_init = linspace(-0.01, 0.01, num_coils); % 高度初始值 (m)

% 优化变量初始化
x0 = [R_init, Z_init];

% 优化目标函数
objective = @(x) -compute_focus_field(x, mu0, I, z_target, r_target, target_radius, num_coils);

% 优化约束（可添加具体限制，如半径范围）
lb = [0.005*ones(1, num_coils), -0.02*ones(1, num_coils)];
ub = [0.1*ones(1, num_coils), 0.02*ones(1, num_coils)];

% 运行优化
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_opt, fval] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);

% 提取优化结果
R_opt = x_opt(1:num_coils);
Z_opt = x_opt(num_coils+1:end);

% 显示优化后的线圈参数
disp('优化后的线圈半径 (m):');
disp(R_opt);
disp('优化后的线圈高度 (m):');
disp(Z_opt);

% 绘制优化后的线圈
figure;
for i = 1:num_coils
    theta = linspace(0, 2*pi, 100);
    x = R_opt(i) * cos(theta);
    y = R_opt(i) * sin(theta);
    z = Z_opt(i) * ones(size(theta));
    plot3(x, y, z, 'LineWidth', 2); hold on;
end
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('优化后的线圈分布');
grid on;
axis equal;

%% 磁场计算函数
function focus_field = compute_focus_field(x, mu0, I, z_target, r_target, target_radius, num_coils)
    R = x(1:num_coils);
    Z = x(num_coils+1:end);

    % 定义中心区域采样点
    [X, Y, Z_target] = meshgrid(linspace(-target_radius, target_radius, 10), ...
                                linspace(-target_radius, target_radius, 10), ...
                                z_target);
    focus_field = 0;
    
    % 毕奥-萨伐尔定律计算磁场
    for k = 1:num_coils
        theta = linspace(0, 2*pi, 100);
        dtheta = theta(2) - theta(1);
        for t = 1:length(theta)
            % 线圈元的位置
            dl = [ -R(k)*sin(theta(t)), R(k)*cos(theta(t)), 0]*dtheta*R(k);
            r_coil = [R(k)*cos(theta(t)), R(k)*sin(theta(t)), Z(k)];
            
            % 对每个采样点计算磁场
            for i = 1:numel(X)
                r = [X(i), Y(i), Z_target(i)];
                r_diff = r - r_coil;
                r_mag = norm(r_diff);
                if r_mag > 1e-6  % 避免奇点
                    dB = mu0 * I * cross(dl, r_diff) / (4*pi*r_mag^3);
                    focus_field = focus_field + norm(dB)^2; % 积分磁场平方
                end
            end
        end
    end
end

% % 网格参数
% Nx = 100; Ny = 100; % 网格大小
% x = linspace(-0.3, 0.3, Nx); % X方向坐标 (m)
% y = linspace(-0.2, 0.2, Ny); % Y方向坐标 (m)
% [XX, YY] = meshgrid(x, y);
% 
% % 初始流函数分布（例如均匀分布）
% psi0 = zeros(Nx, Ny);
% 
% % 目标区域
% target_x = [-0.05, 0.05]; % 中心区域X范围
% target_y = [-0.05, 0.05]; % 中心区域Y范围
% 
% % 优化目标函数
% objective = @(psi) compute_target_field(psi, x, y, target_x, target_y);
% 
% % 优化约束
% lb = -inf * ones(Nx, Ny); % 流函数下界
% ub = inf * ones(Nx, Ny);  % 流函数上界
% 
% % 优化选项
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
% 
% % 执行优化
% [psi_opt, fval] = fmincon(objective, psi0, [], [], [], [], lb, ub, [], options);
% 
% % 绘制优化后的流函数分布
% figure;
% contourf(XX, YY, psi_opt, 20);
% colorbar;
% xlabel('X (m)');
% ylabel('Y (m)');
% title('Optimized Stream Function');
% 
% % 计算目标区域的磁场
% function B = compute_target_field(psi, x, y, target_x, target_y)
%     % 计算电流密度
%     [dpsi_dx, dpsi_dy] = gradient(psi, x, y);
%     Jx = -dpsi_dy; % 电流密度X分量
%     Jy = dpsi_dx;  % 电流密度Y分量
% 
%     % 使用Biot-Savart法则计算磁场
%     mu0 = 4 * pi * 1e-7; % 真空磁导率
%     [XX, YY] = meshgrid(x, y);
%     Bx = mu0 * trapz(trapz(Jx ./ sqrt(XX.^2 + YY.^2)));
%     By = mu0 * trapz(trapz(Jy ./ sqrt(XX.^2 + YY.^2)));
% 
%     % 目标区域磁场积分
%     B = -sum(Bx(:).^2 + By(:).^2); % 最大化磁场强度（取负值）
% end
