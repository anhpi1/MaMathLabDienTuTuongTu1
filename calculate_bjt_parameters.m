%[ZH, AV, ZIN, VCE] = calculate_bjt_parameters(4992450.81, 5304.4050, 13.1758, 5304.4050, 375, 15, 0.675)
function [ZH, AV, ZIN, VCE] = calculate_bjt_parameters(Rb, Rc, Re1, Re2, B, Vcc, Vbe)
    % Tính toán các đại lượng theo tham số truyền vào
    denominator = Rb + (B+1)*(Re1 + Re2);
    Ib = (Vcc - Vbe) / denominator;

    ZH = (26e-3) / ((B+1)*Ib) + Re1;
    AV = -(1/(1/Rc+1/100e3))/ ZH;
    ZIN = 1 / (1/Rb + 1/(B*ZH));
    VCE = Vcc - B*Ib*(Rc + Re1 + Re2);

    % In kết quả kèm ngày hiện tại
    fprintf('A_V test = %.2f\n', AV);
    fprintf('Z_IN test = %.2f Ohm\n', ZIN);
    fprintf('V_CE test= %.2f V\n', VCE);
end
