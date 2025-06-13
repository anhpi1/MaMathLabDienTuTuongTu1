function find_resistor_combo()
    target = 2600;      % Giá trị mục tiêu (Ohm)
    tol = 0.05;            % Sai số cho phép (2%)
    resistors = [1,5,12,22 ,33,39,56,68,75,82,91,100,...
                 220,330,470,560,680,750,820,910,1000,2200,5450,...
                 10000,20000,33000,56000,100000,150000,200000,220000,330000,470000,...
                 1000000,1800000,2200000,3900000,4700000,10000000];  % Ohm

    is_valid = @(R) abs(R - target) <= tol * target;

    combinations = {};

    % --- Hàm tiện ích để thêm tổ hợp ---
    function add_combo(type, R_array, priority)
        val = 0;
        if contains(type, 'parallel')
            val = 1 / sum(1 ./ R_array);
        else
            val = sum(R_array);
        end
        if is_valid(val)
            combinations{end+1} = struct( ...
                'type', type, ...
                'resistors', R_array, ...
                'value', val, ...
                'count', numel(R_array), ...
                'priority', priority ...
            );
        end
    end

    % --------- 1 điện trở đơn ---------
    for i = 1:length(resistors)
        add_combo('single', resistors(i), 1);
    end

    % --------- 2 ~ 4 điện trở song song ---------
    for n = 2:4
        combos = nchoosek(1:length(resistors), n);
        for row = 1:size(combos, 1)
            Rvals = resistors(combos(row, :));
            add_combo(sprintf('%d_parallel', n), Rvals, n);
        end
    end

    % --------- 2 ~ 3 điện trở nối tiếp ---------
    for n = 2:3
        combos = nchoosek(1:length(resistors), n);
        for row = 1:size(combos, 1)
            Rvals = resistors(combos(row, :));
            add_combo(sprintf('%d_series', n), Rvals, n + 3);
        end
    end

    % --------- Chọn tổ hợp tốt nhất ---------
    if isempty(combinations)
        disp('Không tìm được tổ hợp nào thỏa mãn trong sai số cho phép.');
        return;
    end

    best_combo = combinations{1};
    for i = 2:length(combinations)
        c = combinations{i};
        if c.priority < best_combo.priority || ...
           (c.priority == best_combo.priority && c.count < best_combo.count) || ...
           (c.priority == best_combo.priority && c.count == best_combo.count && abs(c.value - target) < abs(best_combo.value - target))
            best_combo = c;
        end
    end

    % --------- Hiển thị kết quả ---------
    fprintf('Tìm được tổ hợp tốt nhất:\n');
    fprintf('- Loại mắc: %s\n', best_combo.type);
    fprintf('- Các điện trở dùng (Ohm): %s\n', mat2str(best_combo.resistors));
    fprintf('- Giá trị đạt được: %.2f Ohm (mục tiêu: %.2f Ohm)\n', best_combo.value, target);
    fprintf('- Số điện trở dùng: %d\n', best_combo.count);
end
