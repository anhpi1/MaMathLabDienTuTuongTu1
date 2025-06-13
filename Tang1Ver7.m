function run_multiple_ga_searches()
    % Constants
    Vcc = 18;
    Vbe = 0.685;
    B = 350;

    % Targets
    Av_target = -7;
    Zin_target = 1e2;
    Zout_target = 1e4;

    %lb = [1e4, 1e3, 0.001, 0.001]; % [Rb, Rc, Re1, Re2]
    %ub = [1e7, 5e3, 500, 1e4];
    % Variable bounds
    lb = [1e4, 1e2, 0.001, 0.001]; % [Rb, Rc, Re1, Re2]
    ub = [1e7, 5e3, 1e4, 1e5];

    % GA options
    options = optimoptions('ga', ...
        'PopulationSize', 500, ...
        'MaxGenerations', 50, ...
        'UseParallel', true, ...
        'HybridFcn', @fmincon, ...
        'FunctionTolerance', 1e-6, ...
        'Display', 'iter');

    % Storage
    best_cost = inf;
    best_solution = [];
    valid_solutions = [];

    attempts = 0;
    runs = 0;
    max_attempts = 5;
    required_runs = 2;
    old_solutions = [];

    while runs < required_runs && attempts < max_attempts
        [x, cost] = ga(@(x)objective(x, Vcc, Vbe, B, Av_target, Zin_target, Zout_target), ...
                       4, [], [], [], [], lb, ub, [], options);

        % Evaluate validity
        [Av, Zin, Zout, Vce] = evaluate_params(x, Vcc, Vbe, B);
        is_valid = Av > Av_target && Zin > Zin_target && Zout < Zout_target && Vce > 8 && Vce < 10;

        % Check duplicates
        is_duplicate = false;
        for i = 1:size(old_solutions,1)
            if all(abs((x - old_solutions(i,:)) ./ (ub - lb)) < 0.01)
                is_duplicate = true;
                break;
            end
        end

        if is_valid && ~is_duplicate
            valid_solutions = [valid_solutions; x, Av, Zin, Zout, Vce]; %#ok<AGROW>
            old_solutions = [old_solutions; x]; %#ok<AGROW>
            runs = runs + 1;
            fprintf('✅ Valid solution found: Run %d\n', runs);
        else
            fprintf('❌ Invalid or duplicate solution. Retrying...\n');
        end

        attempts = attempts + 1;
    end

    % Display best valid solution
    if isempty(valid_solutions)
        fprintf('❌ No valid solution found after %d attempts.\n', attempts);
    else
        % Sort by Av
        [~, idx] = sort(valid_solutions(:,5)); % Sort by Av (smaller is better)
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

function cost = objective(x, Vcc, Vbe, B, Av_target, Zin_target, Zout_target)
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    [Av, Zin, Zout, Vce] = evaluate_params(x, Vcc, Vbe, B);

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
        penalty = penalty + ((Vcc/2 - 0.1*Vcc) - Vce)^2;
    elseif Vce > (Vcc/2) + 0.1*Vcc
        penalty = penalty + (Vce - (Vcc/2 + 0.1*Vcc))^2;
    end

    cost = Av + 1e6 * penalty;
end

function [Av, Zin, Zout, Vce] = evaluate_params(x, Vcc, Vbe, B)
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    I = (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    ZH = 26e-3 / ((B+1)*I) + Re1;
    Av = -Rc / ZH;
    Zin = 1 / (1/Rb + 1/(B*ZH));
    Zout = Rc;
    Vce = Vcc - B * I * (Rc + Re1 + Re2);
end
