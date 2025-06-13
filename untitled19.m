% === Main script ===

% User-defined parameters
Zin_target = 25e3;
Av1_max = -118;
Av2_min = -120;
Vcc = 15;
B = 367;
Zout_max = 10e3;
Vce1_max = 8;
Vce2_min = 7;
Vbe = 0.7; % Typical base-emitter voltage for BJT

% Variable bounds: [Rb, Rc, Re1, Re2]
lb = [4000e3, 1e3, 0.01, 0.001];
ub = [15000e3, 10e3, 10, 10e3];

% Define objective with penalty
objfun = @(x) constraint_penalty(x, Zin_target, Av1_max, Av2_min, ...
                                 Vcc, Vbe, B, Zout_max, Vce1_max, Vce2_min);

% Optimization options
options = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', 500, ...
    'PlotFcn','surrogateoptplot', ...
    'Display','iter');

% Run surrogate optimization
[x, fval] = surrogateopt(objfun, lb, ub, options);

% Extract and display results
Rb_opt = x(1);
Rc_opt = x(2);
Re1_opt = x(3);
Re2_opt = x(4);

fprintf('\nOptimal values found:\n');
fprintf('Rb   = %.2f Ohms\n', Rb_opt);
fprintf('Rc   = %.2f Ohms\n', Rc_opt);
fprintf('Re1  = %.2f Ohms\n', Re1_opt);
fprintf('Re2  = %.2f Ohms\n', Re2_opt);
[ZH, AV, ZIN, VCE] = calculate_bjt_parameters(Rb_opt, Rc_opt, Re1_opt, Re2_opt, 150, 12, 0.7);
function penalty = constraint_penalty(x, Zin_target, Av1_max, Av2_min, ...
                                      Vcc, Vbe, B, Zout_max, Vce1_max, Vce2_min)
    [c, ~] = bjt_constraints(x, Zin_target, Av1_max, Av2_min, ...
                             Vcc, Vbe, B, Zout_max, Vce1_max, Vce2_min);

    % Penalty = sum of squares of violated constraints
    penalty = sum((c .* (c > 0)).^2);
end

function [c, ceq] = bjt_constraints(x, Zin_target, Av1_max, Av2_min, ...
                                    Vcc, Vbe, B, Zout_max, Vce1_max, Vce2_min)
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    % Bias current
    Ib = (Vcc - Vbe) / (Rb + (B + 1) * (Re1 + Re2));

    % Z_H computation
    Z_H = 26e-3 / ((B + 1) * Ib) + Re1;

    % Gain
    Av1 = -Rc / Z_H;
    Av2 = -Rc / Z_H;

    % Input impedance
    Zin_calc = 1 / (1/Rb + 1/(B * Z_H));

    % Collector-emitter voltage
    Vce_common = Vcc - B * Ib * (Rc + Re1 + Re2);

    % Inequality constraints: c(i) <= 0
    c = zeros(7,1);
    c(1) = Av1 - Av1_max;           % Av1 <= Av1_max
    c(2) = Av2_min - Av2;           % Av2 >= Av2_min
    c(3) = Zin_calc - Zin_target;   % Zin <= Zin_target
    c(4) = Vce1_max - Vce_common;   % Vce1 <= max
    c(5) = Vce_common - Vce2_min;   % Vce2 >= min
    c(6) = Rc - Zout_max;           % Rc <= Zout_max
    c(7) = -Z_H;                    % Z_H > 0

    ceq = []; % No equality constraints
end
