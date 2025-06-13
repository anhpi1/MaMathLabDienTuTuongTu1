function transistor_optimization()
    % Constants
    Vcc = 19.5;
    Vbe = 0.685;
    B = 238;

    % Variable bounds
    lb = [10e3, 1e3, 0.001, 0.001];       % Rb, Rc, Re1, Re2 lower bounds
    ub = [10e6, 5e3, 500, 10e3];          % upper bounds

    % GA options
    options = optimoptions('ga', ...
        'PopulationSize',1000, ...
        'MaxGenerations',100, ...
        'UseParallel',true, ...
        'HybridFcn',@fmincon, ...
        'Display','iter', ...
        'PlotFcn',{@gaplotbestf,@gaplotscorediversity});

    % Run GA three times with different random seeds and store results
    results = repmat(struct('x',[],'fval',[],'output',[],'constrviolation',[]),3,1);
    rng('shuffle'); % randomize for first run
    for runIdx=1:3
        % Set different rng seed for each run to get different initial populations
        rng('shuffle');
        [x,fval,exitflag,output] = ga(@(x)objective(x,Vcc,Vbe,B),4,[],[],[],[],lb,ub,[],options);
        % Store results and constraint violation for later comparison
        [~,penalty] = objective(x,Vcc,Vbe,B);
        results(runIdx).x = x;
        results(runIdx).fval = fval;
        results(runIdx).output = output;
        results(runIdx).constrviolation = penalty;
    end

    % Find best feasible solution across runs (penalty == 0 means feasible)
    feasibleIndices = find([results.constrviolation]==0);
    if isempty(feasibleIndices)
        disp('No fully feasible solution found.');
        % Choose best solution with minimum penalty
        [~,bestIdx] = min([results.constrviolation]);
    else
        % Among feasible, choose best (lowest fval = min Av)
        [~,bestIdx] = min([results(feasibleIndices).fval]);
        bestIdx = feasibleIndices(bestIdx);
    end

    bestX = results(bestIdx).x;

    % Calculate final outputs exactly
    [Av, penalty, Zin, Zout, Vce] = compute_params(bestX, Vcc, Vbe, B);

    % Display results
    fprintf('\nBest solution found:\n');
    fprintf('Rb = %.4g Ohms\n', bestX(1));
    fprintf('Rc = %.4g Ohms\n', bestX(2));
    fprintf('Re1= %.4g Ohms\n', bestX(3));
    fprintf('Re2= %.4g Ohms\n', bestX(4));
    fprintf('Av = %.4f\n', Av);
    fprintf('Zin = %.4f Ohms\n', Zin);
    fprintf('Zout = %.4f Ohms\n', Zout);
    fprintf('Vce = %.4f V\n', Vce);
    fprintf('Constraint penalty = %.4f (0 means feasible)\n', penalty);
end

function [obj, penalty] = objective(x, Vcc, Vbe, B)
    % x = [Rb Rc Re1 Re2]
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    % Calculate parameters exactly per given equations
    [Av, ~, Zin, Zout, Vce] = compute_params(x, Vcc, Vbe, B);

    % Objective: minimize Av (most negative is better)
    obj = Av;

    % Penalties for constraint violations (square penalty)
    penalty = 0;

    % Av must be <= -120 (or abs(Av) >= 120)
    if Av > -120
        penalty = penalty + (Av + 120)^2;
    end
    % Zin >= 100k
    if Zin < 100000
        penalty = penalty + (100000 - Zin)^2;
    end
    % Zout <= 10k
    if Zout > 10000
        penalty = penalty + (Zout - 10000)^2;
    end
    % 9 <= Vce <= 10
    if Vce < 9
        penalty = penalty + (9 - Vce)^2;
    elseif Vce > 10
        penalty = penalty + (Vce - 10)^2;
    end

    % Add penalty to objective to force constraints satisfaction
    obj = obj + 1e6 * penalty;
end

function [Av, penalty, Zin, Zout, Vce] = compute_params(x, Vcc, Vbe, B)
    % Unpack variables
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    % Calculate Z_H per given equation:
    % Z_H = 26e-3 / ((B+1)*((Vcc - Vbe)/(Rb + (B+1)*(Re1 + Re2)))) + Re1
    denom = (B+1)*((Vcc - Vbe)/(Rb + (B+1)*(Re1 + Re2)));
    Z_H = 26e-3 / denom + Re1;

    % Av = -Rc / Z_H
    Av = -Rc / Z_H;

    % Zin = parallel combination of Rb and B*Z_H
    Zin = 1 / (1/Rb + 1/(B*Z_H));

    % Vce
    % Vce = Vcc - B * ((Vcc - Vbe)/(Rb + (B+1)*(Re1 + Re2))) * (Rc + Re1 + Re2)
    Vce = Vcc - B * ((Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2))) * (Rc + Re1 + Re2);

    % Zout = Rc (direct)
    Zout = Rc;

    % For convenience
    penalty = 0; % no penalty here, only in objective
end
