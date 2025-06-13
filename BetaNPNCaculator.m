% Các thông số đã cho
Vcc = 12.08;              % Nguồn Vcc (nối emitter) (V)
Vce_measured = 2.76    % Điện áp VCE thực tế cần đạt (V)
Vbe = 0.693;             % Giả sử V_BE (V)
Rb = 500e3;            % Rb = 500k ohm
Rc = 984;              % Rc = 1k ohm
Re = 981;              % Re = 1k ohm

% Hàm phi tuyến cần giải: f(B) = Vce_thuc_te - Vce_tu_cong_thuc = 0
f = @(B) Vcc - B .* ((Vcc - Vbe) ./ (Rb + (B+1)*Re)) .* (Rc + Re) - Vce_measured;

% Đoán ban đầu cho B
B0 = 100;

% Giải phương trình phi tuyến để tìm B
B_solution = fzero(f, B0);

% Hiển thị kết quả
fprintf('Giá trị B tương ứng với Vce = %.2f V là: B = %.2f\n', Vce_measured, B_solution);
