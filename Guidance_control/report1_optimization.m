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

%% plot the contour of the cost function
f = @(x,y) x .* exp(-x.^2 - y.^2) + (x.^2 + y.^2) ./ 20;
x_min = -5; x_max = 5;
y_min = -5; y_max = 5;
n = 100;
[xGrid, yGrid] = meshgrid(linspace(x_min, x_max, n), linspace(y_min, y_max, n));
zGrid = f(xGrid, yGrid);
figure(1);
surfc(xGrid, yGrid, zGrid);
colormap();
colorbar;
xlabel('$x_1$');
ylabel('$x_2$');
title('Contour of the cost function');
axis equal;

figure(2);
contour(xGrid, yGrid, zGrid, 20);
hold on;
colorbar;
xlabel('$x_1$');
ylabel('$x_2$');
title('Contour of the cost function');
axis equal;
hold off;

%% implement the sttepest descent method
gradf = @(x,y) [ (1 - 2*x).*exp(-(x^2 + y^2)) + (1/10)*x;
    -2*x*y*exp(-(x^2 + y^2)) + (1/10)*y ];

H = @(x,y) [ (4*x^2 - 6*x)*exp(-(x^2 + y^2)) + 1/10, (4*(x^2)*y - 2*y)*exp(-(x^2 + y^2));
    (4*(x^2)*y - 2*y)*exp(-(x^2 + y^2)), (4*x*(y^2) - 2*x)*exp(-(x^2 + y^2)) + 1/10 ];