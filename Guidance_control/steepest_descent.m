function [x_hist, f_hist] = steepest_descent(f, grad_f, x0, max_iter, tol)
% STEEPEST_DESCENT_2D   2変数関数に対する最急降下法
%
%  [x_hist, f_hist] = steepest_descent_2d(f, grad_f, x0, max_iter, tol)
%
%  入力:
%    f        : 関数ハンドル (@(x) ...)、x は 2×1 ベクトル
%    grad_f   : 勾配の関数ハンドル (@(x) ...)、2×1 ベクトルを返す
%    x0       : 初期解（2×1 ベクトル）
%    max_iter : 最大反復回数（例: 1000）
%    tol      : 収束判定の閾値（grad ノルム、例: 1e-6）
%
%  出力:
%    x_hist : 各イテレーションの x（2×N マトリクス）
%    f_hist : 各イテレーションの f(x)（1×N ベクトル）

x = x0(:);
x_hist = x;
f_hist = f(x);

for k = 1:max_iter
    g = grad_f(x);
    if norm(g) < tol
        fprintf('Converged steepest descent method at x0 = [%.4f; %.4f]\n', x0(1), x0(2));
        fprintf('  f(x) = %.4f\n', f(x));
        fprintf('  x = [%.4f; %.4f]\n', x(1), x(2));
        fprintf('  k = %d\n', k);
        fprintf('  grad = [%.4f; %.4f]\n', g(1), g(2));
        break;
    end
    if k == max_iter
        fprintf('Do not converged steepest descent method at x0 = [%.4f; %.4f]\n', x0(1), x0(2));
        fprintf('  grad = [%.4f; %.4f]\n', g(1), g(2));
    end

    % バックトラック線形探索（Armijo ルール）
    alpha = 1.0;
    rho   = 0.5;
    c     = 1e-4;
    while f(x - alpha * g) > f(x) - c * alpha * (g' * g)
        alpha = rho * alpha;
    end

    % 更新
    x = x - alpha * g;
    x_hist(:, end+1) = x;
    f_hist(end+1)   = f(x);
end
end
