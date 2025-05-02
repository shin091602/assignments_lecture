close all;  %close all figures
clear;      %clear all variables
clc;        %clear the command terminal
format long
%warning off

% line width
set(0, 'DefaultLineLineWidth', 1.2) % default 0.5pt
set(0, 'DefaultAxesLineWidth', 1.2)
set(0, 'DefaultTextLineWidth', 1.2)

% font size
set(0, 'DefaultTextFontSize', 24)
set(0, 'DefaultAxesFontSize', 24)

% font name
set(0, 'DefaultTextFontName', 'Times New Roman')
set(0, 'DefaultAxesFontName', 'Times New Roman')
set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')

% figure color
set(0, 'DefaultFigureWindowStyle', 'docked');
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardCopy', 'off');

close

% brachistochrone_nd_bvp.m
%--------------------------------------------------
% 最急降下曲線問題を無次元化して bvp4c で解くスクリプト
%--------------------------------------------------


%--- 問題パラメータ ---
g    = 9.81;        % 重力加速度 [m/s^2]
L    = 1.0;         % x-方向の総移動距離 [m]
eps0 = 1e-6;        % 初速度ゼロを避けるための微小パラメータ

%--- メッシュと初期推定値 ---
xmesh     = linspace(0, L, 10000);              
init_guess = [g*eps0; 0];                     % [V; λ] の定数初期推定
solinit   = bvpinit(xmesh, init_guess);

%--- 境界値問題を解く ---
sol = bvp4c(@odefun, @bcfun, solinit);

%--- 結果の取り出しと後処理 ---
x     = sol.x;            % x 座標
V     = sol.y(1, :);      % 速度 V(x)
lam   = sol.y(2, :);      % ラグランジュ乗数 λ(x)
theta = asin(g * lam);    % 最適角度 θ(x) : sinθ = g·λ
y     = cumtrapz(x, tan(theta));  % y(x) = ∫ tanθ dx

%--- プロット ---
figure;
plot(x, -y, 'LineWidth', 2);
grid on;
xlabel('x [m]');
ylabel('y [m]');
title('Brachistochrone via bvp4c');

function dYdx = odefun(x, Y)
% Y = [V; λ]
g   = 9.81;
V   = Y(1);
lam = Y(2);
th  = asin(g * lam);
dVdx   = (g / V) * tan(th);
dlamdx = (1 / V^2) * (1 / cos(th) + g * lam * tan(th));
dYdx   = [dVdx; dlamdx];
end

function res = bcfun(Y0, YL)
% 境界条件: V(0)=g·eps0,  λ(L)=0
g    = 9.81;
eps0 = 1e-6;
res = [ Y0(1) - g*eps0;
    YL(2)         ];
end
