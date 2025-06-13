% Parameters
Vcc = 19.5;
Vbe = 0.7;
B = 238;
VT = 26e-3; % Thermal voltage

% Lower and upper bounds
lb = [1e4, 500, 1, 500];      % [Rb, Rc, Re1, Re2]
ub = [1e7, 10000, 200, 2e3];

% Objective function (minimize Rc)
obj_fun = @(x) x(2); % x = [Rb Rc Re1 Re2]

% Nonlinear constraints
nonlcon = @(x) constraints(x, Vcc, Vbe, B, VT);

% Run genetic algorithm
opts = optimoptions('ga','Display','iter','MaxGenerations',300);
[x_opt, fval] = ga(obj_fun, 4, [], [], [], [], lb, ub, nonlcon, opts);

% Extract optimal values
Rb = x_opt(1);
Rc = x_opt(2);
Re1 = x_opt(3);
Re2 = x_opt(4);

% Recalculate for output
Ic = B * ((Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2)));
Zin = 1 / (1/Rb + 1/(B*(VT / ((B+1)*((Vcc - Vbe)/(Rb + (B+1)*(Re1 + Re2)))) + Re1)));
Av = Rc / (VT / ((B+1)*((Vcc - Vbe)/(Rb + (B+1)*(Re1 + Re2)))) + Re1);
Vce = Vcc - B*((Vcc - Vbe)/(Rb + (B+1)*(Re1 + Re2))) * (Rc + Re1 + Re2);

% Display
fprintf('\nOptimal Solution:\n');
fprintf('Rb = %.2f Ohms\n', Rb);
fprintf('Rc = %.2f Ohms\n', Rc);
fprintf('Re1 = %.2f Ohms\n', Re1);
fprintf('Re2 = %.2f Ohms\n', Re2);
fprintf('Av = %.4f\n', Av);
fprintf('Vce = %.4f V\n', Vce);
fprintf('Zin = %.2f Ohms\n', Zin);
fprintf('Ic = %.4f A\n', Ic);
[ZH, AV, ZIN, VCE] = calculate_bjt_parameters(Rb, Rc, Re1, Re2, B, Vcc, Vbe);

function [c, ceq] = constraints(x, Vcc, Vbe, B, VT)
    Rb = x(1); Rc = x(2); Re1 = x(3); Re2 = x(4);

    Ib = (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Ic = B * Ib;
    
    Zin = 1 / (1/Rb + 1/(B*(VT / ((B+1)*Ib) + Re1)));
    Av = Rc / (VT / ((B+1)*Ib) + Re1);
    Vce = Vcc - Ic * (Rc + Re1 + Re2);

    % Inequality constraints (c <= 0)
    c = [...
        100000 - Zin;      % Zin > 10000
        9 - Vce;          % Vce > 9
        Vce - 10;         % Vce < 10
        58 - Av;          % Av > 14
        Av - 60;          % Av < 16
        Ic - 0.4          % Ic < 0.4
    ];

    ceq = []; % No equality constraints
end


