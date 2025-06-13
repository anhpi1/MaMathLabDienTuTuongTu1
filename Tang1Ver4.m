function ga_bjt_optimization
    % Tham số cố định
    Vcc = 19.5;
    Vbe = 0.685;
    B = 238;

    % Biên biến số (giới hạn tìm kiếm)
    lb = [1e3, 1e2, 1e-3, 1e-3];       % [Rb, Rc, Re1, Re2]
    ub = [1e8, 5e3, 500, 1e4];

    opts = optimoptions('ga',...
        'MaxGenerations', 100,...
        'PopulationSize', 300,...
        'UseParallel', true,...
        'Display','iter',...
        'PlotFcn', {@gaplotbestf});

    foundValid = false;
    while ~foundValid
        [x,fval,exitflag,output] = ga(@(params) objFun(params, Vcc, Vbe, B),4,[],[],[],[],lb,ub,[],opts);

        % Tính lại các giá trị tại nghiệm tìm được
        [Av, Zin, Zout, Vce] = calculateAll(x, Vcc, Vbe, B);

        % Kiểm tra điều kiện ràng buộc
        if (Av >= -120) && (Zin >= 1e5) && (Zout <= 1e4) && (Vce >= 9) && (Vce <= 10)
            foundValid = true;
            fprintf('Found valid solution:\n');
            fprintf('Rb = %.4g Ohm\n', x(1));
            fprintf('Rc = %.4g Ohm\n', x(2));
            fprintf('Re1 = %.4g Ohm\n', x(3));
            fprintf('Re2 = %.4g Ohm\n', x(4));
            fprintf('Av = %.4f\n', Av);
            fprintf('Zin = %.4g Ohm\n', Zin);
            fprintf('Zout = %.4g Ohm\n', Zout);
            fprintf('Vce = %.4f V\n', Vce);
            
        else
            fprintf('No valid solution this run, restarting GA...\n');
        end
    end
end

function cost = objFun(params, Vcc, Vbe, B)
    % params = [Rb, Rc, Re1, Re2]
    Rb = params(1);
    Rc = params(2);
    Re1 = params(3);
    Re2 = params(4);

    % Tính Z_H theo công thức nguyên bản
    denom = (B+1)*(Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Z_H = (26e-3) / denom + Re1;

    % Tính Av
    Av = -Rc / Z_H;

    % Tính Zin
    Zin = 1 / (1/Rb + 1/(B * Z_H));

    % Tính Vce
    current_term = B * (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Vce = Vcc - current_term * (Rc + Re1 + Re2);

    % Tính Zout = R_E1 (cách hiểu phổ biến, không rõ bạn muốn gì hơn)
    Zout = Re1;

    % Biến mục tiêu: chọn Av (nhỏ nhất)
    % Cộng penalty lớn nếu vi phạm ràng buộc
    penalty = 0;
    if Av < -120
        penalty = penalty + 1e6 * abs(Av + 120);
    end
    if Zin < 1e5
        penalty = penalty + 1e6 * (1e5 - Zin);
    end
    if Zout > 1e4
        penalty = penalty + 1e6 * (Zout - 1e4);
    end
    if Vce < 9
        penalty = penalty + 1e6 * (9 - Vce);
    end
    if Vce > 10
        penalty = penalty + 1e6 * (Vce - 10);
    end

    cost = Av + penalty;  % tối thiểu Av với penalty nếu vi phạm ràng buộc
end

function [Av, Zin, Zout, Vce] = calculateAll(params, Vcc, Vbe, B)
    Rb = params(1);
    Rc = params(2);
    Re1 = params(3);
    Re2 = params(4);

    denom = (B+1)*(Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Z_H = (26e-3) / denom + Re1;
    Av = -Rc / Z_H;
    Zin = 1 / (1/Rb + 1/(B * Z_H));
    current_term = B * (Vcc - Vbe) / (Rb + (B+1)*(Re1 + Re2));
    Vce = Vcc - current_term * (Rc + Re1 + Re2);
    Zout = Re1; % như đề bài không rõ Zout là gì, lấy Re1 làm Zout

end
