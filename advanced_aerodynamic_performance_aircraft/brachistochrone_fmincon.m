% brachistochrone_fmincon.m
%--------------------------------------------------
% Brachistochrone を fmincon で直接転写解法
% 固定端点 (0,0) -> (L,H)
%--------------------------------------------------
close all; clear; clc;

%% 問題パラメータ
g = 9.81;    % 重力加速度 [m/s^2]
L = 1.0;     % 水平方向距離 [m]
H = 0.5;     % 垂直落下量 [m] （始点 y=0, 終点 y=H>0）
N = 100;     % 分割数（適宜増減可）

%% グリッド生成
x  = linspace(0, L, N)';      % x_i
dx = x(2) - x(1);             % Δx

%% 初期推定：直線
y0 = linspace(0, H, N)';      % y_i の初期推定

%% 目的関数：時間の積分を離散化
timefun = @(y) ...
    sum( sqrt(1 + ((diff(y)/dx).^2)) .* dx ./ ...
    sqrt( g * (y(1:end-1) + y(2:end)) ) );

%% 線形等式制約：始点・終点固定
Aeq = zeros(2, N);
Aeq(1,1) = 1;    beq(1) = 0;   % y(1)=0
Aeq(2,N) = 1;    beq(2) = H;   % y(N)=H

%% 線形不等式制約：単調増加 y(i+1)-y(i)>=0
M = N-1;
A = zeros(M, N);
for i = 1:M
    A(i,i)   = -1;
    A(i,i+1) =  1;
end
b = zeros(M,1);

%% 境界：y >= 0
lb = zeros(N,1);
ub = [];

%% fmincon オプション
opts = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxFunEvals',1e5, ...
    'MaxIter',1000);

%% 最適化実行
[y_opt, T_opt] = fmincon( timefun, y0, A, b, Aeq, beq, lb, ub, [], opts );

%% 結果プロット
figure;
plot(x, y_opt, '-o','LineWidth',1.5);
grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Brachistochrone via fmincon (T = %.4f s)', T_opt));
