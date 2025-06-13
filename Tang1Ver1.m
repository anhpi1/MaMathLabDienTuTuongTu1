% BJT Amplifier Design Optimization using Genetic Algorithm
% Goal: Maximize Av but Av <= 120
% Constraints: Zin >= 100e3, Zout <= 10e3, 9 <= Vce <= 10
% Given parameters:
Vcc = 19.5;
Vbe = 0.685;
beta = 238;

% Variable bounds:
lb = [10e3, 1e3, 1e-3, 1e-3];      % [Rb, Rc, Re1, Re2]
ub = [1e6, 10e3, 500, 10e3];     % [Rb, Rc, Re1, Re2]

% Objective: maximize Av, so minimize -Av
objFun = @(x) -computeAv(x, Vcc, Vbe, beta);

% Nonlinear constraints function
nonlcon = @(x) bjtConstraints(x, Vcc, Vbe, beta);

% GA options
options = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 200, ...
    'UseParallel', false);

solutionFound = false;
while ~solutionFound
    [x_opt, fval, exitflag] = ga(objFun, 4, [], [], [], [], lb, ub, nonlcon, options);
    if exitflag > 0
        solutionFound = true;
    else
        disp('No valid solution found, retrying GA...');
    end
end

% Extract optimized values
Rb_opt  = x_opt(1);
Rc_opt  = x_opt(2);
Re1_opt = x_opt(3);
Re2_opt = x_opt(4);

% Compute final metrics
Zin_opt  = computeZin(x_opt, Vcc, Vbe, beta);
Zout_opt = Rc_opt;                    % Approximate output impedance
Av_opt   = computeAv(x_opt, Vcc, Vbe, beta);
Vce_opt  = computeVce(x_opt, Vcc, Vbe, beta);

% Display results
fprintf('Optimized Design:\n');
fprintf('Rb = %.2f kΩ\n', Rb_opt/1e3);
fprintf('Rc = %.2f kΩ\n', Rc_opt/1e3);
fprintf('Re1 = %.2f Ω\n', Re1_opt);
fprintf('Re2 = %.2f Ω\n', Re2_opt);
fprintf('Zin = %.2f kΩ\n', Zin_opt/1e3);
fprintf('Zout = %.2f kΩ\n', Zout_opt/1e3);
fprintf('Av = %.2f\n', Av_opt);
fprintf('Vce = %.2f V\n', Vce_opt);
[ZH, AV, ZIN, VCE] = calculate_bjt_parameters(Rb_opt, Rc_opt, Re1_opt, Re2_opt, B, Vcc, Vbe);

%% Function definitions
function Av = computeAv(x, Vcc, Vbe, beta)
    Rb  = x(1);
    Rc  = x(2);
    Re1 = x(3);
    Re2 = x(4);
    Zh = (26e-3) / ((beta+1)*( (Vcc-Vbe)/(Rb + (beta+1)*(Re1+Re2)) )) + Re1;
    Av = -Rc / Zh;
end

function Zin = computeZin(x, Vcc, Vbe, beta)
    Rb  = x(1);
    Re1 = x(3);
    Re2 = x(4);
    Zh = (26e-3) / ((beta+1)*( (Vcc-Vbe)/(Rb + (beta+1)*(Re1+Re2)) )) + Re1;
    Zin = (1/Rb + 1/(beta * Zh))^(-1);
end

function Vce = computeVce(x, Vcc, Vbe, beta)
    Rb  = x(1);
    Rc  = x(2);
    Re1 = x(3);
    Re2 = x(4);
    Ib = (Vcc - Vbe) / (Rb + (beta+1)*(Re1+Re2));
    Ic = beta * Ib;
    Vce = Vcc - Ic * (Rc + Re1 + Re2);
end

function [c, ceq] = bjtConstraints(x, Vcc, Vbe, beta)
    % Inequality constraints c(x) <= 0
    Zin  = computeZin(x, Vcc, Vbe, beta);
    Zout = x(2);           % Rc approx
    Av   = computeAv(x, Vcc, Vbe, beta);
    Vce  = computeVce(x, Vcc, Vbe, beta);
    c = [100e3 - Zin; ...    % Zin >= 100e3 -> 100e3 - Zin <=0
         Zout - 10e3; ...    % Zout <= 10e3 -> Zout - 10e3 <=0
         9 - Vce; ...        % Vce >= 9    -> 9 - Vce <=0
         Vce - 10; ...       % Vce <= 10   -> Vce - 10 <=0
         Av - 120];          % Av <= 120   -> Av -120 <=0
    ceq = [];
end
