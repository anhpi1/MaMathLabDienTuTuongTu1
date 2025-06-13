% Parameters
Vcc = 18;
Vbe = 1.228;
B = 23973;

objFun = @(x) computeZout(x(1), x(2), Vcc, Vbe, B);

nonlcon = @(x) constraints(x(1), x(2), Vcc, Vbe, B);

lb = [5e5, 0.0001];
ub = [1e6, 0.01];

options = optimoptions('ga','Display','iter','MaxGenerations',100);

[x_opt,zout_opt] = ga(objFun, 2, [], [], [], [], lb, ub, nonlcon, options);

Rb = x_opt(1);
Re = x_opt(2);
Ic = (Vcc - Vbe) / (Rb + (B+1)*Re);
Vce = Vcc - B * Ic * Re;
Zin = computeZin(Rb, Re, Vcc, Vbe, B);
Zout = zout_opt;

fprintf('Optimal Rb = %.3f Ohm\n', Rb);
fprintf('Optimal Re = %.6f Ohm\n', Re);
fprintf('Vce = %.4f V\n', Vce);
fprintf('Zin = %.2f Ohm\n', Zin);
fprintf('Zout = %.4f Ohm\n', Zout);

% === Functions ===
function Zin = computeZin(Rb, Re, Vcc, Vbe, B)
    Ic = (Vcc - Vbe) / (Rb + (B+1)*Re);
    term = (26e-3 / (B+1)) * Ic + Re;
    Zin = 1/Rb + B * (1 / term);
end

function Zout = computeZout(Rb, Re, Vcc, Vbe, B)
    numerator = (B+1)*(Vcc - Vbe);
    denominator = Rb + (B+1)*Re;
    Zout_inv = (1/Re) + (numerator / denominator) / (26e-3);
    Zout = 1 / Zout_inv;
end

function [c, ceq] = constraints(Rb, Re, Vcc, Vbe, B)
    Ic = (Vcc - Vbe) / (Rb + (B+1)*Re);
    Vce = Vcc - B * Ic * Re;
    Zin = computeZin(Rb, Re, Vcc, Vbe, B);
    c = [...
        19.45 - Vce;      % Vce > 9
        Vce - 19.5;     % Vce < 10
        20000 - Zin   % Zin > 20000
    ];
    ceq = [];
end
