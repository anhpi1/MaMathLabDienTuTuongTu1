function optimize_bjt_amplifier()
% ------------------------------------------------------------
%  Genetic-Algorithm optimisation of a BJT common-emitter stage
% ------------------------------------------------------------
%  Variables (decision vector x) : [R1  R2  RC  RE]  (Ohm)
% ------------------------------------------------------------

%% ----------- User-defined constants ------------------------
Vcc      = 15;        % [V]   Supply
Vbe      = 0.675;     % [V]   Base-emitter drop
beta     = 275;       %       DC current gain
RL       = 10e3;      % [Ohm] External load
VT       = 25e-3;     % [V]   Thermal voltage (≈25 mV @ 25 °C)

%% ----------- Search bounds ---------------------------------
%            R1       R2        RC      RE
lb = [   1e3,    1e3,     100,     100];        % lower
ub = [   1e6,    1e6,   50e3,    5e3];          % upper

%% ----------- GA options ------------------------------------
opts = optimoptions('ga', ...
    'HybridFcn',        @fmincon, ...
    'PopulationSize',   150, ...
    'MaxGenerations',   200, ...
    'FunctionTolerance',1e-6, ...
    'ConstraintTolerance',1e-6, ...
    'Display',          'iter');

%% ----------- Objective & constraints -----------------------
obj     = @(x) objectiveFun (x, Vcc, Vbe, beta, RL, VT);
nonlcon = @(x) nonlinearCon(x, Vcc, Vbe, beta, RL, VT);

%% ----------- Run GA ----------------------------------------
[x_opt, fval] = ga(obj, 4, [], [], [], [], lb, ub, nonlcon, opts);

%% ----------- Report result ---------------------------------
[R1,R2,RC,RE] = deal(x_opt(1),x_opt(2),x_opt(3),x_opt(4));
fprintf('\n===== Optimised Resistor Values =====\n');
fprintf('R1 = %.2f  Ohm\n', R1);
fprintf('R2 = %.2f  Ohm\n', R2);
fprintf('RC = %.2f  Ohm\n', RC);
fprintf('RE = %.2f  Ohm\n', RE);
fprintf('Max |Av| found = %.2f\n', -fval);

verifySolution(R1,R2,RC,RE,Vcc,Vbe,beta,RL,VT);
end
% =====================================================================
% Objective: minimise  -|Av|   (i.e. maximise |Av|)
% =====================================================================
function cost = objectiveFun(x,Vcc,Vbe,beta,RL,VT)
    [~,~,~,~,Av] = smallSignal(x,Vcc,Vbe,beta,RL,VT);
    cost = -abs(Av);
end
% =====================================================================
% Non-linear inequality constraints   c(x) ≤ 0
% =====================================================================
function [c,ceq] = nonlinearCon(x,Vcc,Vbe,beta,RL,VT)
    [Ic,Vce,Zin,~,Av] = smallSignal(x,Vcc,Vbe,beta,RL,VT);

    c = [ ...
        7.0 - Vce      ;      % Vce ≥ 7
        Vce - 8.0      ;      % Vce ≤ 8
        1e-3 - Ic      ;      % Ic ≥ 1 mA
        10e3 - Zin     ;      % Zin ≥ 10 kΩ
        abs(Av) - 50   ;      % |Av| ≤ 50
        40 - abs(Av)          % |Av| ≥ 40
    ];
    ceq = [];                 % no equality constraints
end
% =====================================================================
% Helper: compute bias & small-signal parameters for a given x
% =====================================================================
function [Ic,Vce,Zin,re,Av] = smallSignal(x,Vcc,Vbe,beta,RL,VT)
    R1 = x(1);  R2 = x(2);  RC = x(3);  RE = x(4);

    % Thevenin at the base
    Rth  = parallel(R1,R2);
    Eth  = R2*Vcc/(R1+R2);

    % DC currents/voltages
    Ib   = (Eth - Vbe) / (Rth + (beta+1)*RE);
    Ic   = beta * Ib;
    Vce  = Vcc - Ic*(RC + RE);

    % Small-signal
    re   = VT / Ic;
    Av   = -parallel(RL,RC) / re;
    Zin  = parallel(parallel(R1,R2), beta*re);   % CORRECT formula
end
% =====================================================================
% Verify final solution & assert constraints
% =====================================================================
function verifySolution(R1,R2,RC,RE,Vcc,Vbe,beta,RL,VT)
    [Ic,Vce,Zin,re,Av] = smallSignal([R1,R2,RC,RE],Vcc,Vbe,beta,RL,VT);

    fprintf('\n===== Verified Circuit Parameters =====\n');
    fprintf('IC  = %.4f  A\n', Ic);
    fprintf('VCE = %.2f   V\n', Vce);
    fprintf('re  = %.2f  Ohm\n', re);
    fprintf('Zin = %.2f  kOhm\n', Zin/1e3);
    fprintf('Zout= %.2f  Ohm\n', RC);
    fprintf('|Av| = %.2f\n', abs(Av));

    % Hard assertions (fail loudly if violated)
    assert(Vce >= 7 && Vce <= 8,   'VCE constraint violated');
    assert(Ic  >= 1e-3,            'IC  constraint violated');
    assert(Zin >= 10e3,            'Zin constraint violated');
    assert(abs(Av) >= 40 && abs(Av) <= 50, '|Av| constraint violated');
end
% =====================================================================
% Utility: two-resistor parallel
% =====================================================================
function rp = parallel(Ra,Rb)
    rp = (Ra.*Rb)./(Ra+Rb);
end
