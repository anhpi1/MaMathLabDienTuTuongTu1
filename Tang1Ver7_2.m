function run_multiple_ga_searches()
    % Constants
    Vcc = 15;
    Vbe = 0.685;
    B = 375;

    % Targets
    Av_target = -240;
    Zin_target = 15e3;
    Zout_target = 1e4;
    %lb = [1e4, 1e2, 0.001, 5.6e3]; % [Rb, Rc, Re1, Re2]
    %ub = [2e7, 20e3, 1e7, 1e7];
    % Variable bounds [Rb, Rc, Re1, Re2]

    lb = [9900e3, 9.9e3, 0.001, 1]; % [Rb, Rc, Re1, Re2]
    ub = [11000e3, 10.1e3, 1e4, 1e7];

    % GA options
    options = optimoptions('ga', ...
        'PopulationSize', 100000, ...
        'MaxGenerations', 100, ...
        'UseParallel', true, ...
        'HybridFcn', @fmincon, ...
        'FunctionTolerance', 1e-6, ...
        'Display', 'iter');

    % Storage
    valid_solutions = [];  % Mảng lưu các nghiệm hợp lệ: [x, Av, Zin, Zout, Vce]
    all_solutions = [];    % Mảng lưu tất cả nghiệm (bao gồm cả penalty)
    old_valids = [];       % Lưu lại các x của các nghiệm hợp lệ đã đc chọn

    attempts = 0;
    runs = 0;
    max_attempts = 1;     % Số lần tối đa lặp lại thuật toán
    required_runs = 10;     % Số lần tối thiểu lặp lại thuật toán

    while runs < required_runs && attempts < max_attempts
        attempts = attempts + 1;

        [x, cost] = ga(@(x) objective(x, Vcc, Vbe, B, Av_target, Zin_target, Zout_target), ...
                       4, [], [], [], [], lb, ub, [], options);

        % Evaluate candidate x (từng nghiệm tính được)
        [Av, Zin, Zout, Vce] = evaluate_params(x, Vcc, Vbe, B);
        penalty = compute_penalty(Av, Zin, Zout, Vce, Av_target, Zin_target, Zout_target, Vcc);
        is_valid = (penalty == 0);

        % --- Kiểm tra trùng lặp với các cá thể đã có (ngoại trừ 200 cá thể tốt nhất) ---
        % Nếu có đủ nghiệm hợp lệ, tách ra top 200 (theo thứ tự ưu tiên; vì nghiệm hợp lệ có penalty=0 nên ta sắp xếp theo Av, ví dụ)
        if ~isempty(valid_solutions)
            num_best = min(20, size(valid_solutions, 1));
            % Sắp xếp theo giá trị Av (ở cột 5) – nhỏ hơn (trường hợp âm) được ưu tiên
            [~, idx_sorted] = sort(valid_solutions(:,5));  
            best_valid = valid_solutions(idx_sorted(1:num_best), 1:4);  % chỉ lấy phần thiết kế x

            % Lấy ra danh sách các nghiệm hợp lệ (x) không thuộc top 200
            all_valid = valid_solutions(:, 1:4);
            others = [];
            for j = 1:size(all_valid, 1)
                if ~is_member_approx(all_valid(j,:), best_valid, lb, ub, 0.01)
                    others = [others; all_valid(j,:)];  %#ok<AGROW>
                end
            end
        else
            others = [];
        end

        % Kiểm tra nếu x trùng với bất kỳ cá thể trong "others" (theo bán kính 1% của khoảng tìm kiếm)
        is_duplicate = false;
        if ~isempty(others)
            for i = 1:size(others, 1)
                if all( abs((x - others(i,:)) ./ (ub - lb)) < 0.01 )
                    is_duplicate = true;
                    break;
                end
            end
        end

        % Lưu toàn bộ nghiệm (bao gồm cả không hợp lệ)
        all_solutions = [all_solutions; x, Av, Zin, Zout, Vce, penalty];   %#ok<AGROW>

        % Nếu nghiệm hợp lệ và không bị trùng (x không nằm trong những cá thể ngoài top 200)
        if is_valid && ~is_duplicate
            valid_solutions = [valid_solutions; x, Av, Zin, Zout, Vce];   %#ok<AGROW>
            old_valids = [old_valids; x];   %#ok<AGROW>
            runs = runs + 1;
            fprintf('✅ Valid solution found: Run %d\n', runs);
        else
            fprintf('❌ Invalid or duplicate solution. Retrying...\n');
        end
    end

    if isempty(valid_solutions)
        fprintf('❌ No valid solution found after %d attempts.\n', attempts);
        % In nghiệm không hợp lệ với penalty nhỏ nhất
        [~, idx] = min(all_solutions(:, end)); % chọn nghiệm có penalty nhỏ nhất
        best_invalid = all_solutions(idx, :);
        fprintf('\n⚠️ Closest Invalid Solution (minimal penalty):\n');
        fprintf('Rb = %.2f\n', best_invalid(1));
        fprintf('Rc = %.2f\n', best_invalid(2));
        fprintf('Re1 = %.4f\n', best_invalid(3));
        fprintf('Re2 = %.4f\n', best_invalid(4));
        fprintf('Av = %.2f\n', best_invalid(5));
        fprintf('Zin = %.2f\n', best_invalid(6));
        fprintf('Zout = %.2f\n', best_invalid(7));
        fprintf('Vce = %.2f\n', best_invalid(8));
        fprintf('Penalty = %.2f\n', best_invalid(9));
    else
        % Sắp xếp nghiệm hợp lệ theo Av (ưu tiên giá trị càng âm tốt hơn)
        [~, idx] = sort(valid_solutions(:,5));
        best_solution = valid_solutions(idx(1), :);
        fprintf('\n🎯 Best Valid Solution:\n');
        fprintf('Rb = %.2f\n', best_solution(1));
        fprintf('Rc = %.2f\n', best_solution(2));
        fprintf('Re1 = %.4f\n', best_solution(3));
        fprintf('Re2 = %.4f\n', best_solution(4));
        fprintf('Av = %.2f\n', best_solution(5));
        fprintf('Zin = %.2f\n', best_solution(6));
        fprintf('Zout = %.2f\n', best_solution(7));
        fprintf('Vce = %.2f\n', best_solution(8));
    end
end

% Hàm so sánh xem một vector candidate có "gần" với một trong các vector có trong mảng group hay không
function isMember = is_member_approx(candidate, group, lb, ub, tol)
    % tol được tính theo tỉ lệ của (ub-lb)
    isMember = false;
    for k = 1:size(group,1)
        if all( abs((candidate - group(k,:)) ./ (ub - lb)) < tol )
            isMember = true;
            return;
        end
    end
end

function cost = objective(x, Vcc, Vbe, B, Av_target, Zin_target, Zout_target)
    [Av, Zin, Zout, Vce] = evaluate_params(x, Vcc, Vbe, B);
    penalty = compute_penalty(Av, Zin, Zout, Vce, Av_target, Zin_target, Zout_target, Vcc);
    cost = Av + 1e6 * penalty;
end

function penalty = compute_penalty(Av, Zin, Zout, Vce, Av_target, Zin_target, Zout_target, Vcc)
    penalty = 0;
    if Av < Av_target
        penalty = penalty + (Av - Av_target)^2;
    end
    if Zin < Zin_target
        penalty = penalty + (Zin_target - Zin)^2;
    end
    if Zout > Zout_target
        penalty = penalty + (Zout - Zout_target)^2;
    end
    if Vce < (Vcc/2) - 0.1*Vcc
        penalty = penalty + (((Vcc/2)-0.1*Vcc) - Vce)^2;
    elseif Vce > (Vcc/2) + 0.1*Vcc
        penalty = penalty + (Vce - ((Vcc/2)+0.1*Vcc))^2;
    end
end

function [Av, Zin, Zout, Vce] = evaluate_params(x, Vcc, Vbe, B)
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    I = (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    ZH = 26e-3 / ((B+1)*I) + Re1;
    Av = -(1/(1/Rc+1/100e3)) / ZH;
    Zin = 1 / (1/Rb + 1/(B*ZH));
    Zout = Rc;
    Vce = Vcc - B * I * (Rc + Re1 + Re2);
end
