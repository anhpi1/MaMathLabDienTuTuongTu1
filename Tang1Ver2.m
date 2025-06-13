% transistor_ga_parfor.m

% Problem constants
Vcc = 19.5;
Vbe = 1.4;
B = 90000;
kT = 26e-3;  % Thermal voltage

% Boundaries for Rb, Rc, Re1, Re2
lb = [1, 1,   0.001,  0.001];
ub = [1e8,  1e8,  1e8,    1e8];

% GA Options
opts = optimoptions("ga", ...
    "PopulationSize", 300, ...
    "MaxGenerations", 100, ...
    "MaxStallGenerations", 100, ...
    "InitialPopulationRange", [lb; ub], ...
    "Display", "iter", ...
    "UseParallel", true);

% Wrapper to restart GA until solution is found
solutionFound = false;
while ~solutionFound
    [x, fval, exitflag, output] = ga(@(x) penalized_gain(x, Vcc, Vbe, B, kT), ...
        4, [], [], [], [], lb, ub, [], opts);

    [Av_abs, Zin, Zout, Vce] = metrics(x, Vcc, Vbe, B, kT);

    if Av_abs <= 120 && Zin > 1e5 && Zout < 1e4 && Vce > 9 && Vce < 10
        solutionFound = true;
        fprintf("Found valid solution:\n");
        fprintf("Rb  = %.2f Ohm\n", x(1));
        fprintf("Rc  = %.2f Ohm\n", x(2));
        fprintf("Re1 = %.6f Ohm\n", x(3));
        fprintf("Re2 = %.6f Ohm\n", x(4));
        fprintf("Zin  = %.2f Ohm\n", Zin);
        fprintf("Zout = %.2f Ohm\n", Zout);
        fprintf("|Av| = %.2f\n", Av_abs);
        fprintf("Vce  = %.2f V\n", Vce);
        [ZH, AV, ZIN, VCE] = calculate_bjt_parameters(x(1), x(2), x(3), x(4), B, Vcc, Vbe);
    else
        fprintf("GA did not converge to feasible point – restarting search...\n\n");
    end
end

% Objective with constraint penalty
function f = penalized_gain(x, Vcc, Vbe, B, kT)
    [Av_abs, Zin, Zout, Vce] = metrics(x, Vcc, Vbe, B, kT);
    penalty = 0;
    if Av_abs > 120;      penalty = penalty + 1e5; end
    if Zin < 1e5;         penalty = penalty + 1e5; end
    if Zout > 1e4;        penalty = penalty + 1e5; end
    if Vce < 9 || Vce > 10; penalty = penalty + 1e5; end
    f = -Av_abs + penalty;
end

% Calculate all dependent variables
function [Av_abs, Zin, Zout, Vce] = metrics(x, Vcc, Vbe, B, kT)
    Rb  = x(1);
    Rc  = x(2);
    Re1 = x(3);
    Re2 = x(4);

    Ib = (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Zh = kT / ((B+1)*Ib) + Re1;
    Av_abs = abs(-Rc / Zh);
    Zin = 1 / (1/Rb + 1/(B*Zh));
    Zout = Rc;  % Simple approximation assuming CE stage
    Ic = B * Ib;
    Vce = Vcc - Ic * (Rc + Re1 + Re2);
end