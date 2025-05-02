% DynamicProgrammingForBrachistochroneSimplified.m
% 動的計画法でブラキストクローン問題を解き、その曲線 y(x) をプロット

clear; close all; clc;

%% 問題設定
g       = 9.81;      % 重力加速度 [m/s^2]
x_end   = 2.0;       % 水平方向移動距離 [m]
y_end   = 1.2;       % 垂直落差 [m]
nx      = 100;        % x 方向分割数
ny      = 800;        % y 方向分割数
delta_x = x_end / (nx - 1);
delta_y = y_end / (ny - 1);
% 動的計画法で最適経路を求める
[y_opt, x_grid] = BrachistochroneDP(delta_x, delta_y, nx, ny);
% 結果プロット
figure;
y_plot = y_end - y_opt;
plot(x_grid, y_opt, 'b-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse');   % ← ここで y 軸を反転
xlabel('x [m]');
ylabel('y(x) [m]');
title('Brachistochrone Curve via Dynamic Programming');
grid on;
% 関数定義：動的計画法ソルバ
function [y_output, x_grid] = BrachistochroneDP(dx, dy, nx, ny)
g = 9.81;
% グリッド
x_grid = (0:nx-1) * dx;
y_grid = (0:ny-1) * dy;
% コストと前向きインデックスの初期化
f_all   = inf(nx, ny);
prev_idx = zeros(nx, ny);

% スタートから x=dx への初期ステップ
for j = 2:ny
    v_avg = sqrt(g/2) * (sqrt(y_grid(j)) + 0);
    f_all(2,j) = sqrt(dx^2 + y_grid(j)^2) / max(v_avg, realmin);
    prev_idx(2,j) = 1;
end

% メインの動的計画法ループ
for i = 3:nx
    for j = 2:ny
        best = inf; best_k = 1;
        for k = 1:ny
            v_avg = sqrt(g/2) * (sqrt(y_grid(j)) + sqrt(y_grid(k)));
            t_step = sqrt(dx^2 + (y_grid(j)-y_grid(k))^2) / max(v_avg, realmin);
            cost = f_all(i-1,k) + t_step;
            if cost < best
                best = cost;
                best_k = k;
            end
        end
        f_all(i,j) = best;
        prev_idx(i,j) = best_k;
    end
end

% 経路のバックトラック
y_output = zeros(1, nx);
idx = ny;                % 終端は y= y_end のインデックス
y_output(end) = y_grid(idx);
for i = nx:-1:2
    idx = prev_idx(i, idx);
    y_output(i-1) = y_grid(idx);
end
end
