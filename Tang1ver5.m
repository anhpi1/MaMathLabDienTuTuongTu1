function optimize_bjt_circuit
    % Tham số cố định
    Vcc = 15;
    Vbe = 0.685;
    B = 350;
    Av_target= -120;      % Av phải >= -120 (tức Av không được âm quá lớn)
    Zout_target= 1e4;     % Zout < 10k
    Zin_target=1e5;       % Zin > 100k

    % Biên cho biến quyết định
    % Rb(10k 10M), Rc(1k 5k), Re1(0.001 500), Re2(0.001 10k)
    lb = [1e6, 1e3, 1e-3, 1e-3];
    ub = [10e6, 5e3, 500, 1e4];

    % Options GA với hybrid function (fmincon local search), song song
    options = optimoptions('ga', ...
        'PopulationSize', 1000, ...      % giảm kích thước để nhanh hơn
        'MaxGenerations', 100, ...       % tăng số thế hệ để đủ thời gian chạy
        'UseParallel', true, ...
        'HybridFcn', @fmincon, ...
        'Display', 'iter');

    solution_found = false;
    while ~solution_found
        [x,fval,exitflag,output] = ga(@(x)objective(x, Vcc, Vbe, B, Av_target, Zin_target, Zout_target), 4, ...
            [], [], [], [], lb, ub, [], options);

        % Tính lại các đại lượng với nghiệm tìm được
        [Av, Zin, Zout, Vce] = circuit_params(x, Vcc, Vbe, B);

        % Kiểm tra điều kiện ràng buộc
        if (Av >= Av_target) && (Zin > Zin_target) && (Zout < Zout_target) && ...
           (Vce >= (Vcc/2)-0.1*Vcc) && (Vce <= (Vcc/2)+0.1*Vcc)
            solution_found = true;
            % Hiển thị kết quả
            fprintf('\nSolution found:\n');
            fprintf('Rb = %.3g Ohm\n', x(1));
            fprintf('Rc = %.3g Ohm\n', x(2));
            fprintf('Re1= %.3g Ohm\n', x(3));
            fprintf('Re2= %.3g Ohm\n', x(4));
            fprintf('Av = %.4f\n', Av);
            fprintf('Zin = %.3g Ohm\n', Zin);
            fprintf('Zout= %.3g Ohm\n', Zout);
            fprintf('Vce = %.4f V\n', Vce);
            [ZH, AV, ZIN, VCE] = calculate_bjt_parameters(x(1), x(2), x(3), x(4), B, Vcc, Vbe);
        else
            fprintf('\nSolution not found, retrying...\n');
        end
    end
end

function cost = objective(x, Vcc, Vbe, B, Av_target, Zin_target, Zout_target)
    % Tính các đại lượng
    [Av, Zin, Zout, Vce] = circuit_params(x, Vcc, Vbe, B);

    % Penalty square function cho các ràng buộc
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
    if Vce < (Vcc/2)-0.1*Vcc
        penalty = penalty + (((Vcc/2)-0.1*Vcc) - Vce)^2;
    elseif Vce > (Vcc/2)+0.1*Vcc
        penalty = penalty + (Vce - ((Vcc/2)+0.1*Vcc))^2;
    end
    cost = Av + 1e6*penalty; % penalty nhân hằng số lớn
end

function [Av, Zin, Zout, Vce] = circuit_params(x, Vcc, Vbe, B)
    Rb = x(1);
    Rc = x(2);
    Re1 = x(3);
    Re2 = x(4);

    % Tính Z_H theo công thức
    denominator = (B+1)*(Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Z_H = 26e-3 / denominator + Re1;

    % Tính Av
    Av = -Rc / Z_H;

    % Tính Zin
    Zin = 1 / (1/Rb + 1/(B*Z_H));

    % Tính Vce
    I_factor = B * (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Vce = Vcc - I_factor * (Rc + Re1 + Re2);

    % Tính Zout
    Zout = Rc;
end
