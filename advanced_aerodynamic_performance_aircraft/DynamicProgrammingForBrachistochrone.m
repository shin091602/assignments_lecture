% DynamicProgrammingForBrachistochrone.m
% Corrected by ChatGPT: fixes for string handling, initial prev_output, Octave/MATLAB compatibility, and legend handling

clear variables;
close all;
clc;

% Load Octave packages if running in Octave
if exist('OCTAVE_VERSION','builtin')
  pkg load control;
  pkg load signal;
end

format long e;

% Variant selection (use single quotes for char arrays)
variant = 'y(2.0)=1.2';
%variant = 'y(1.0)=1.0';

% Define end coordinates based on variant
if strcmp(variant, 'y(1.0)=1.0')
  x_end = 1.0;
  y_end = 1.0;
elseif strcmp(variant, 'y(2.0)=1.2')
  x_end = 2.0;
  y_end = 1.2;
else
  error('Unknown variant: %s', variant);
end

% Discretization parameters
x_point_nmb = 26; % number of x grid points
y_point_nmb = 51; % number of y grid points
Y_end_index = y_point_nmb;

% Run dynamic programming solver
[y_output, f_output, x_grid, y_grid] = BrachistochroneCurveWithDynamicProgramming( ...
  x_end/(x_point_nmb-1), y_end/(y_point_nmb-1), ...
  x_point_nmb, y_point_nmb, Y_end_index);

% Optimal time from DP
time_dynamic_programming_in_seconds_opt = f_output(end);

% Known optimal time from variation method (for comparison)
if strcmp(variant, 'y(1.0)=1.0')
  time_variation_in_seconds_opt = 0.58288;
elseif strcmp(variant, 'y(2.0)=1.2')
  time_variation_in_seconds_opt = 0.8006;
end
percentMistake = abs(time_dynamic_programming_in_seconds_opt - time_variation_in_seconds_opt) * 100 / time_variation_in_seconds_opt;


% Plot DP results
figure;
plot(x_grid, y_output, 'b-', x_grid, f_output, 'r-');
legend({'y output', 'f output'});
xlabel('x [m]');
ylabel('y / time [m or s]');
grid on;

% Variation Method support curves
theta = 0:0.01:2*pi;
figure;
f1 = 1 - cos(theta);
f2 = (y_end / x_end) .* (theta - sin(theta));
plot(theta, f1, theta, f2);
title('Brachistochrone Variation Method Support');
legend({'1 - cos(	heta)', '(Y_{end}/X_{end})(	heta - sin(	heta))'});
xlabel('	heta [rad]');
ylabel('Function value');
grid on;

% Compute r_end via iteration or graphical value
if strcmp(variant, 'y(1.0)=1.0')
  r_end = IterationFunction(3, 25);
elseif strcmp(variant, 'y(2.0)=1.2')
  r_end = 3.234;
end

const1 = 2 * y_end / (1 - cos(r_end));
const2 = 2 * x_end / (r_end - sin(r_end));
r = 0:0.1:r_end;
x_var = (const2/2) * (r - sin(r));
y_var = (const1/2) * (1 - cos(r));
x_var(end+1) = x_end;
y_var(end+1) = y_end;

% Plot variation curve
figure;
plot(x_var, y_var);
legend(sprintf('Variation: y(0)=0, y(%.1f)=%.1f', x_end, y_end));
grid on;
axis ij;
xt = x_grid(1:5:end);
xticks = xt;
set(gca, 'XTick', x_grid(1:5:end), 'XAxisLocation', 'top');
xlabel('x [m]');;
ylabel('y [m]');

% Combined plot
figure;
plot(x_var, y_var, 'b-', x_grid, y_output, 'r--');
legend({'Variation', 'Dynamic Programming'});
grid on;
axis ij;
set(gca, 'XTick', x_grid(1:5:end), 'XAxisLocation', 'top');
xlabel('x [m]');
ylabel('y [m]');

%-------------------------------------------------------------------------%
% Function: BrachistochroneCurveWithDynamicProgramming
function [y_output, f_output, x_grid, y_grid] = BrachistochroneCurveWithDynamicProgramming(delta_x, delta_y, x_point_nmb, y_point_nmb, Y_end_index)
g = 9.81;
x_grid = (0:x_point_nmb-1) * delta_x;
y_grid = (0:y_point_nmb-1) * delta_y;
f_all = inf(x_point_nmb, y_point_nmb);
prev_idx = zeros(x_point_nmb, y_point_nmb);

% Initial step from start (x=1)
for j = 2:y_point_nmb
  v_avg = sqrt(g/2) * (sqrt(y_grid(j)) + sqrt(0));
  if v_avg > realmin
    f_all(2,j) = sqrt(delta_x^2 + (y_grid(j))^2) / v_avg;
  end
  prev_idx(2,j) = 1;
end

% Main DP loops
for i = 3:(x_point_nmb-1)
  for j = 2:y_point_nmb
    f_curr = inf;
    for k = 1:y_point_nmb
      v_avg = sqrt(g/2) * (sqrt(y_grid(j)) + sqrt(y_grid(k)));
      if v_avg > realmin
        t_step = sqrt(delta_x^2 + (y_grid(j)-y_grid(k))^2) / v_avg;
        f_temp = f_all(i-1,k) + t_step;
      else
        f_temp = inf;
      end
      if f_temp < f_curr
        f_curr = f_temp;
        prev_idx(i,j) = k;
      end
    end
    f_all(i,j) = f_curr;
  end
end

% Last step to end
f_curr = inf;
for k = 2:y_point_nmb
  v_avg = sqrt(g/2) * (sqrt(y_grid(Y_end_index)) + sqrt(y_grid(k)));
  if v_avg > realmin
    t_step = sqrt(delta_x^2 + (y_grid(Y_end_index)-y_grid(k))^2) / v_avg;
    f_temp = f_all(x_point_nmb-1,k) + t_step;
  else
    f_temp = inf;
  end
  if f_temp < f_curr
    f_curr = f_temp;
    prev_idx(x_point_nmb, Y_end_index) = k;
  end
end
f_all(x_point_nmb, Y_end_index) = f_curr;

% Backtrack path
y_output = zeros(1, x_point_nmb);
f_output = zeros(1, x_point_nmb);
idx = Y_end_index;
y_output(end) = y_grid(idx);
f_output(end) = f_all(end, idx);
for i = (x_point_nmb-1):-1:2
  idx = prev_idx(i+1, idx);
  y_output(i) = y_grid(idx);
  f_output(i) = f_all(i, idx);
end
y_output(1) = 0;
f_output(1) = 0;
end

%-------------------------------------------------------------------------%
% Function: IterationFunction
function x_end = IterationFunction(x_start, iteration_nmb)
x_prev = x_start;
for i = 1:iteration_nmb
  x_prev = 1 - cos(x_prev) + sin(x_prev);
end
x_end = x_prev;
end
