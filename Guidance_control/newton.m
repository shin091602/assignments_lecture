function [x_hist, f_hist, eig_hist] = newton(f, grad_f, hess_f, x0, max_iter, tol)
% NEWTON_WITH_EIG   ニュートン法＋ヘッセ固有値履歴出力
%
%  [x_hist, f_hist, eig_hist] = newton_with_eig(...)
%
%  入力:
%    f, grad_f, hess_f, x0, max_iter, tol は従来の newton_2d と同じ
%  出力:
%    x_hist   : 各イテレーションの x（2×N マトリクス）
%    f_hist   : 各イテレーションの f(x)（1×N ベクトル）
%    eig_hist : 各イテレーションの H の固有値（2×N マトリクス）

% 初期化
x        = x0(:);
x_hist   = x;
f_hist   = f(x);
eig_hist = zeros(2,1);

% 初期点での固有値
H0       = hess_f(x);
eig_hist(:,1) = eig(H0);

for k = 1:max_iter
    g = grad_f(x);
    if norm(g) < tol
        fprintf('Converged Newton method at x0 = [%.4f; %.4f]\n', x0(1), x0(2));
        fprintf('  f(x) = %.4f\n', f(x));
        fprintf('  x = [%.4f; %.4f]\n', x(1), x(2));
        fprintf('  k = %d\n', k);
        fprintf('  grad = [%.4f; %.4f]\n', g(1), g(2));
        break;
    end
    if k == max_iter
        fprintf('Do not converged Newton method at x0 = [%.4f; %.4f]\n', x0(1), x0(2));
        fprintf('  grad = [%.4f; %.4f]\n', g(1), g(2));
    end

    Hmat = hess_f(x);
    % 固有値を求めて保存
    D = eig(Hmat);
    eig_hist(:, end+1) = D;

    % d = - Hmat \ g;
    if all(D > 0)
        d = - Hmat \ g;
    else
        d = - g;
    end

    % 更新
    x = x + d;
    x_hist(:, end+1) = x;
    f_hist(end+1)    = f(x);
end
end
