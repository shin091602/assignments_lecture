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

%% gradient vector and hessian matrix and initial points
f = @(x) x(1).*exp(-x(1).^2 - x(2).^2) + ( x(1).^2 + x(2).^2 )/20;

gradf = @(x) [ ...
    (1 - 2*x(1).^2) .* exp(-(x(1).^2 + x(2).^2)) + (1/10) * x(1);
    -2 * x(1) .* x(2) .* exp(-(x(1).^2 + x(2).^2)) + (1/10) * x(2)
    ];

H = @(x) [ ...
    (4*x(1).^3 - 6*x(1)) .* exp(-(x(1).^2 + x(2).^2)) + 1/10,   (4*x(1).^2 .* x(2) - 2*x(2)) .* exp(-(x(1).^2 + x(2).^2));
    (4*x(1).^2 .* x(2) - 2*x(2)) .* exp(-(x(1).^2 + x(2).^2)),  (4*x(1) .* x(2).^2 - 2*x(1)) .* exp(-(x(1).^2 + x(2).^2)) + 1/10
    ];
x0_1     = [-0.4; 0.6];
x0_2     = [1.0; 0.1];
x0_3     = [1.6; 1.9];
%% implement the gradient descent method
tol = 1e-4;
[X_1, F] = steepest_descent(f, gradf, x0_1, 100, tol);
figure(3);
plot(F);
xlabel('Iteration');
ylabel('Cost function value');
title('gradient descent method x0 = [-0.4; 0.6]');
grid on

[X_2, F] = steepest_descent(f, gradf, x0_2, 200, tol);
figure(4);
plot(F);
xlabel('Iteration');
ylabel('Cost function value');
title('gradient descent method (x0 = [1.0; 0.1])');
grid on

[X_3, F] = steepest_descent(f, gradf, x0_3, 100, tol);
figure(5);
plot(F);
xlabel('Iteration');
ylabel('Cost function value');
title('gradient descent method (x0 = [1.6; 1.9])');
grid on

%% 軌跡を x–y 平面上にプロット
figure(6);
hold on;
plot(X_1(1, :), X_1(2, :), '-o','DisplayName', 'x0 = [-0.4; 0.6]', 'MarkerFaceColor','auto');
plot(X_2(1, :), X_2(2, :), '-s','DisplayName', 'x0 = [1.0; 0.1]', 'MarkerFaceColor','auto');
plot(X_3(1, :), X_3(2, :), '-^','DisplayName', 'x0 = [1.6; 1.9]','MarkerFaceColor','auto');
contour(xGrid, yGrid, zGrid,20,'handlevisibility','off');
plot(x0_1(1), x0_1(2), 'ko', 'MarkerSize',10, 'LineWidth',2,'HandleVisibility','off');
plot(x0_2(1), x0_2(2), 'ks', 'MarkerSize',10, 'LineWidth',2,'HandleVisibility','off');
plot(x0_3(1), x0_3(2), 'k^', 'MarkerSize',10, 'LineWidth',2,'HandleVisibility','off');
xlabel('$x_1$');
ylabel('$x_2$');
xlim([-2.1, 2.1]);
ylim([-2.1, 2.1]);
title('Gradient Descent Trajectories in x-y Plane');
legend('Location','best');
grid on;
hold off;

%% implement the newton method
[X_1, F, eh] = newton(f, gradf, H, x0_1, 100, tol);
figure(7);
plot(F);
xlabel('Iteration');
ylabel('Cost function value');
title('Newton method (x0 = [-0.4; 0.6])');
grid on

figure(10);
hold on;
plot(eh(1, :), '-o','LineWidth', 1.5, 'DisplayName', '$\lambda_1$');
plot(eh(2, :), '-s','LineWidth', 1.5, 'DisplayName', '$\lambda_2$');
xlabel('Iteration');
ylabel('Eigenvalue');
title('Eigenvalues of Hessian (x0 = [-0.4; 0.6])');
legend('Location','best');
grid on;
hold off;

[X_2, F, eh] = newton(f, gradf, H, x0_2, 400, tol);
figure(8);
plot(F);
xlabel('Iteration');
ylabel('Cost function value');
title('Newton method (x0 = [1.0; 0.1])');
grid on

figure(11);
hold on;
plot(eh(1, :), '-o','LineWidth', 1.5, 'DisplayName', '$\lambda_1$');
plot(eh(2, :), '-s','LineWidth', 1.5, 'DisplayName', '$\lambda_2$');
xlabel('Iteration');
ylabel('Eigenvalue');
title('Eigenvalues of Hessian (x0 = [1.0; 0.1])');
legend('Location','best');
grid on;
hold off;

[X_3, F, eh] = newton(f, gradf, H, x0_3, 400, tol);
figure(9);
plot(F);
xlabel('Iteration');
ylabel('Cost function value');
title('Newton method (x0 = [1.6; 1.9])');
grid on

figure(12);
hold on;
plot(eh(1, :), '-o','LineWidth', 1.5, 'DisplayName', '$\lambda_1$');
plot(eh(2, :), '-s','LineWidth', 1.5, 'DisplayName', '$\lambda_2$');
xlabel('Iteration');
ylabel('Eigenvalue');
title('Eigenvalues of Hessian (x0 = [1.6; 1.9])');
legend('Location','best');
grid on;
hold off;

%% 軌跡を x–y 平面上にプロット
figure(13);
hold on;
plot(X_1(1, :), X_1(2, :), '-o','DisplayName', 'x0 = [-0.4; 0.6]', 'MarkerFaceColor','auto');
plot(X_2(1, :), X_2(2, :), '-s','DisplayName', 'x0 = [1.0; 0.1]', 'MarkerFaceColor','auto');
plot(X_3(1, :), X_3(2, :), '-^','DisplayName', 'x0 = [1.6; 1.9]','MarkerFaceColor','auto');
contour(xGrid, yGrid, zGrid,20,'handlevisibility','off');
plot(x0_1(1), x0_1(2), 'ko', 'MarkerSize',10, 'LineWidth',2,'handlevisibility','off');
plot(x0_2(1), x0_2(2), 'ks', 'MarkerSize',10, 'LineWidth',2,'HandleVisibility','off');
plot(x0_3(1), x0_3(2), 'k^', 'MarkerSize',10, 'LineWidth',2,'HandleVisibility','off');
xlabel('$x_1$');
ylabel('$x_2$');
xlim([-3, 3]);
ylim([-3, 3]);
title('Newton method Trajectories in x-y Plane');
legend('Location','best');
grid on;
hold off;





%% implement the fmincon method
nonlcon = @(x) deal( ...
    x(1)*x(2)/2 ...
    + (x(1)+2).^2 ...
    + (x(2)-2).^2/2 ...
    - 2, ...   % c(x) <= 0
    [] ...    % ceq(x) = []
    );
x0 = [-0.4; 0.6];

% Interior‐Point
opts_IP = optimoptions('fmincon', ...
    'Algorithm',               'interior-point', ...
    'Display',                 'iter-detailed', ...
    'MaxIterations',            100, ...
    'OptimalityTolerance',    1e-8, ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn',              @(x,lambda) H(x) ...
    );

[x_ip, fval_ip, exitflag_ip, output_ip, lambda_ip] = fmincon(@(x) deal(f(x), gradf(x)), x0, [], [], [], [], [], [], nonlcon, opts_IP);

%  Active‐Set
opts_AS = optimoptions('fmincon',...
    'Algorithm',            'active-set',...
    'Display',              'iter',...
    'MaxFunctionEvaluations',1e4);
[x_as,f_as] = fmincon(f, x0, [],[],[],[],[],[], nonlcon, opts_AS);

%% plot the result
figure(14);
hold on;
h1 = plot(x_ip(1), x_ip(2), 'go', ...
    'MarkerSize',10, 'LineWidth',2, ...
    'DisplayName','Interior-Point Solution');

% Active-Set 解
h2 = plot(x_as(1), x_as(2), 'b*', ...
    'MarkerSize',10, 'LineWidth',2, ...
    'DisplayName','Active-Set Solution');

% コスト関数の等高線
h3 = contour(xGrid, yGrid, zGrid, 20, ...
    'handlevisibility','off');

% 制約境界 g=0 の線
gGrid = xGrid.*yGrid/2 ...
    + (xGrid+2).^2 ...
    + (yGrid-2).^2/2 ...
    - 2;
h4 = contour(xGrid, yGrid, gGrid, [0 0], ...
    'LineColor','k', 'LineWidth',2, ...
    'DisplayName','Constraint Boundary');

% 凡例を表示
legend('Location','best');
grid on;
xlabel('$x_1$'); ylabel('$x_2$');
title('Optimization Results with Constraint');


