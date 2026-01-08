global yesgenealogy;
yesgenealogy = 1;  % 1 - строить генеалогию, 0 - нет

% Базовые параметры
distribution_s = 'const';
r = 0;      % Частота рекомбинации
M = 3;      % Количество точек рекомбинации
s0 = 0.1;   % Базовая сила отбора
L = 300;    % Длина генома (количество локусов)
N = 500;   % Размер популяции
tf = 150;   % Время моделирования
f0 = 0;     % Начальная частота благоприятных аллелей
muL = 0.01; % Общая частота мутаций
run = 1;    % номер запуска

%[TMRCA, adapt_data] = recomb_train_fig1(distribution_s, r, M, s0, L, N, tf, f0, muL, run, 'test');
%[TMRCA, adapt_data] = recomb_train(distribution_s, r, M, s0, L, N, tf, f0, muL, run, 'test');
%[TMRCA, adapt_data] = recomb_train_original(distribution_s, r, M, s0, L, N, tf, f0, muL, run);


% Создаем структуру для хранения результатов
results = struct();

%% Серия 1: Влияние размера популяции (N)
disp('=== СЕРИЯ 1: Влияние размера популяции (N) ===');
N_values = [1000, 5000]; % Размеры популяции из НИР

for i = 1:length(N_values)
    current_N = N_values(i);
    [TMRCA, adapt_data] = recomb_train(distribution_s, r, M, s0, L, current_N, tf, f0, muL, run, sprintf('exp_N_%d', i));
    
    % Сохраняем результаты
    results.N(i).N = current_N;
    results.N(i).TMRCA = TMRCA;
    results.N(i).adapt_data = adapt_data;
    
    % Извлекаем данные
    if isstruct(adapt_data) && isfield(adapt_data, 'V_an')
        results.N(i).V_an = adapt_data.V_an;
    else
        % Если структура не содержит V_an, вычисляем аналитическую скорость
        results.N(i).V_an = compute_analytical_velocity(current_N, s0, L, f0, muL);
    end
    
    % Извлекаем k_av
    if isstruct(adapt_data) && isfield(adapt_data, 'k_av')
        k_av_data = adapt_data.k_av;
    elseif isstruct(adapt_data)
        field_names = fieldnames(adapt_data);
        for j = 1:length(field_names)
            field_data = adapt_data.(field_names{j});
            if isnumeric(field_data) && length(field_data) > 1
                k_av_data = field_data;
                break;
            end
        end
    elseif isnumeric(adapt_data)
        k_av_data = adapt_data;
    else
        k_av_data = [];
    end
    
    % Сохраняем данные о k_av
    if ~isempty(k_av_data)
        results.N(i).k_av = k_av_data;
        results.N(i).k_final = k_av_data(end);
        results.N(i).f_final = k_av_data(end) / L;
        
        % Вычисляем численную скорость адаптации
        if length(k_av_data) > 1
            time_steps = 0:length(k_av_data)-1; % Время в поколениях
            % Используем линейную регрессию на последней трети данных
            t_start = round(2*length(time_steps)/3);
            if length(time_steps(t_start:end)) > 2
                p = polyfit(time_steps(t_start:end), k_av_data(t_start:end), 1);
                results.N(i).V_num = p(1);
            else
                % Если мало точек, используем среднюю скорость
                results.N(i).V_num = (k_av_data(end) - k_av_data(1)) / (time_steps(end) - time_steps(1));
            end
        else
            results.N(i).V_num = NaN;
        end
    else
        results.N(i).k_final = NaN;
        results.N(i).f_final = NaN;
        results.N(i).V_num = NaN;
    end
    
    % Выводим результаты
    fprintf('N = %d:\n', current_N);
    fprintf('  TMRCA = %.2f\n', TMRCA);
    fprintf('  Численная скорость адаптации V_num = %.4f\n', results.N(i).V_num);
    fprintf('  Аналитическая скорость адаптации V_an = %.4f\n', results.N(i).V_an);
    if ~isnan(results.N(i).k_final)
        fprintf('  Финальное число благоприятных аллелей k_кон = %.2f\n', results.N(i).k_final);
        fprintf('  Финальная частота благоприятных аллелей f_кон = %.4f\n\n', results.N(i).f_final);
    end
end

%% Серия 2: Влияние числа локусов (L)
disp('=== СЕРИЯ 2: Влияние числа локусов (L) ===');
L_values = [100, 200, 300]; % Число локусов из НИР

for i = 1:length(L_values)
    current_L = L_values(i);
    [TMRCA, adapt_data] = recomb_train(distribution_s, r, M, s0, current_L, N, tf, f0, muL, run, sprintf('exp_L_%d', i));
    
    % Сохраняем результаты
    results.L(i).L = current_L;
    results.L(i).TMRCA = TMRCA;
    results.L(i).adapt_data = adapt_data;
    
    % Вычисляем аналитическую скорость
    results.L(i).V_an = compute_analytical_velocity(N, s0, current_L, f0, muL);
    
    % Извлекаем k_av
    if isstruct(adapt_data) && isfield(adapt_data, 'k_av')
        k_av_data = adapt_data.k_av;
    elseif isstruct(adapt_data)
        field_names = fieldnames(adapt_data);
        for j = 1:length(field_names)
            field_data = adapt_data.(field_names{j});
            if isnumeric(field_data) && length(field_data) > 1
                k_av_data = field_data;
                break;
            end
        end
    elseif isnumeric(adapt_data)
        k_av_data = adapt_data;
    else
        k_av_data = [];
    end
    
    % Сохраняем данные о k_av
    if ~isempty(k_av_data)
        results.L(i).k_av = k_av_data;
        results.L(i).k_final = k_av_data(end);
        results.L(i).f_final = k_av_data(end) / current_L;
        
        % Вычисляем численную скорость адаптации
        if length(k_av_data) > 1
            time_steps = 0:length(k_av_data)-1;
            t_start = round(2*length(time_steps)/3);
            if length(time_steps(t_start:end)) > 2
                p = polyfit(time_steps(t_start:end), k_av_data(t_start:end), 1);
                results.L(i).V_num = p(1);
            else
                results.L(i).V_num = (k_av_data(end) - k_av_data(1)) / (time_steps(end) - time_steps(1));
            end
        else
            results.L(i).V_num = NaN;
        end
    else
        results.L(i).k_final = NaN;
        results.L(i).f_final = NaN;
        results.L(i).V_num = NaN;
    end
    
    % Выводим результаты
    fprintf('L = %d:\n', current_L);
    fprintf('  TMRCA = %.2f\n', TMRCA);
    fprintf('  Численная скорость адаптации V_num = %.4f\n', results.L(i).V_num);
    fprintf('  Аналитическая скорость адаптации V_an = %.4f\n', results.L(i).V_an);
    if ~isnan(results.L(i).k_final)
        fprintf('  Финальное число благоприятных аллелей k_кон = %.2f\n', results.L(i).k_final);
        fprintf('  Финальная частота благоприятных аллелей f_кон = %.4f\n\n', results.L(i).f_final);
    end
end

%% Серия 3: Влияние частоты мутаций (µ_L)
disp('=== СЕРИЯ 3: Влияние частоты мутаций (µ_L) ===');
muL_values = [0.05, 0.1, 0.2]; % Частоты мутаций из НИР

for i = 1:length(muL_values)
    current_muL = muL_values(i);
    [TMRCA, adapt_data] = recomb_train(distribution_s, r, M, s0, L, N, tf, f0, current_muL, run, sprintf('exp_muL_%d', i));
    
    % Сохраняем результаты
    results.muL(i).muL = current_muL;
    results.muL(i).TMRCA = TMRCA;
    results.muL(i).adapt_data = adapt_data;
    
    % Вычисляем аналитическую скорость
    results.muL(i).V_an = compute_analytical_velocity(N, s0, L, f0, current_muL);
    
    % Извлекаем k_av
    if isstruct(adapt_data) && isfield(adapt_data, 'k_av')
        k_av_data = adapt_data.k_av;
    elseif isstruct(adapt_data)
        field_names = fieldnames(adapt_data);
        for j = 1:length(field_names)
            field_data = adapt_data.(field_names{j});
            if isnumeric(field_data) && length(field_data) > 1
                k_av_data = field_data;
                break;
            end
        end
    elseif isnumeric(adapt_data)
        k_av_data = adapt_data;
    else
        k_av_data = [];
    end
    
    % Сохраняем данные о k_av
    if ~isempty(k_av_data)
        results.muL(i).k_av = k_av_data;
        results.muL(i).k_final = k_av_data(end);
        results.muL(i).f_final = k_av_data(end) / L;
        
        % Вычисляем численную скорость адаптации
        if length(k_av_data) > 1
            time_steps = 0:length(k_av_data)-1;
            t_start = round(2*length(time_steps)/3);
            if length(time_steps(t_start:end)) > 2
                p = polyfit(time_steps(t_start:end), k_av_data(t_start:end), 1);
                results.muL(i).V_num = p(1);
            else
                results.muL(i).V_num = (k_av_data(end) - k_av_data(1)) / (time_steps(end) - time_steps(1));
            end
        else
            results.muL(i).V_num = NaN;
        end
    else
        results.muL(i).k_final = NaN;
        results.muL(i).f_final = NaN;
        results.muL(i).V_num = NaN;
    end
    
    % Выводим результаты
    fprintf('µ_L = %.2f:\n', current_muL);
    fprintf('  TMRCA = %.2f\n', TMRCA);
    fprintf('  Численная скорость адаптации V_num = %.4f\n', results.muL(i).V_num);
    fprintf('  Аналитическая скорость адаптации V_an = %.4f\n', results.muL(i).V_an);
    if ~isnan(results.muL(i).k_final)
        fprintf('  Финальное число благоприятных аллелей k_кон = %.2f\n', results.muL(i).k_final);
        fprintf('  Финальная частота благоприятных аллелей f_кон = %.4f\n\n', results.muL(i).f_final);
    end
end

%% Серия 4: Влияние силы отбора (s0)
disp('=== СЕРИЯ 4: Влияние силы отбора (s0) ===');
s0_values = [0.05, 0.1, 0.2]; % Меньшая, базовая и большая сила отбора

for i = 1:length(s0_values)
    current_s0 = s0_values(i);
    [TMRCA, adapt_data] = recomb_train(distribution_s, r, M, current_s0, L, N, tf, f0, muL, run, sprintf('exp_s0_%d', i));
    
    % Сохраняем результаты
    results.s0(i).s0 = current_s0;
    results.s0(i).TMRCA = TMRCA;
    results.s0(i).adapt_data = adapt_data;
    
    % Вычисляем аналитическую скорость
    results.s0(i).V_an = compute_analytical_velocity(N, current_s0, L, f0, muL);
    
    % Извлекаем k_av
    if isstruct(adapt_data) && isfield(adapt_data, 'k_av')
        k_av_data = adapt_data.k_av;
    elseif isstruct(adapt_data)
        field_names = fieldnames(adapt_data);
        for j = 1:length(field_names)
            field_data = adapt_data.(field_names{j});
            if isnumeric(field_data) && length(field_data) > 1
                k_av_data = field_data;
                break;
            end
        end
    elseif isnumeric(adapt_data)
        k_av_data = adapt_data;
    else
        k_av_data = [];
    end
    
    % Сохраняем данные о k_av
    if ~isempty(k_av_data)
        results.s0(i).k_av = k_av_data;
        results.s0(i).k_final = k_av_data(end);
        results.s0(i).f_final = k_av_data(end) / L;
        
        % Вычисляем численную скорость адаптации
        if length(k_av_data) > 1
            time_steps = 0:length(k_av_data)-1;
            t_start = round(2*length(time_steps)/3);
            if length(time_steps(t_start:end)) > 2
                p = polyfit(time_steps(t_start:end), k_av_data(t_start:end), 1);
                results.s0(i).V_num = p(1);
            else
                results.s0(i).V_num = (k_av_data(end) - k_av_data(1)) / (time_steps(end) - time_steps(1));
            end
        else
            results.s0(i).V_num = NaN;
        end
    else
        results.s0(i).k_final = NaN;
        results.s0(i).f_final = NaN;
        results.s0(i).V_num = NaN;
    end
    
    % Выводим результаты
    fprintf('s0 = %.2f:\n', current_s0);
    fprintf('  TMRCA = %.2f\n', TMRCA);
    fprintf('  Численная скорость адаптации V_num = %.4f\n', results.s0(i).V_num);
    fprintf('  Аналитическая скорость адаптации V_an = %.4f\n', results.s0(i).V_an);
    if ~isnan(results.s0(i).k_final)
        fprintf('  Финальное число благоприятных аллелей k_кон = %.2f\n', results.s0(i).k_final);
        fprintf('  Финальная частота благоприятных аллелей f_кон = %.4f\n\n', results.s0(i).f_final);
    end
end

%% Сводная таблица результатов
disp('=== СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ ===');
fprintf('\nСерия 1: Влияние размера популяции (N)\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', 'N', 'TMRCA', 'V_num', 'V_an', 'k_кон', 'f_кон');
for i = 1:length(N_values)
    fprintf('%-10d %-10.2f %-10.4f %-10.4f %-10.2f %-10.4f\n', ...
        results.N(i).N, results.N(i).TMRCA, results.N(i).V_num, results.N(i).V_an, ...
        results.N(i).k_final, results.N(i).f_final);
end

fprintf('\nСерия 2: Влияние числа локусов (L)\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', 'L', 'TMRCA', 'V_num', 'V_an', 'k_кон', 'f_кон');
for i = 1:length(L_values)
    fprintf('%-10d %-10.2f %-10.4f %-10.4f %-10.2f %-10.4f\n', ...
        results.L(i).L, results.L(i).TMRCA, results.L(i).V_num, results.L(i).V_an, ...
        results.L(i).k_final, results.L(i).f_final);
end

fprintf('\nСерия 3: Влияние частоты мутаций (µ_L)\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', 'µ_L', 'TMRCA', 'V_num', 'V_an', 'k_кон', 'f_кон');
for i = 1:length(muL_values)
    fprintf('%-10.2f %-10.2f %-10.4f %-10.4f %-10.2f %-10.4f\n', ...
        results.muL(i).muL, results.muL(i).TMRCA, results.muL(i).V_num, results.muL(i).V_an, ...
        results.muL(i).k_final, results.muL(i).f_final);
end

fprintf('\nСерия 4: Влияние силы отбора (s0)\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', 's0', 'TMRCA', 'V_num', 'V_an', 'k_кон', 'f_кон');
for i = 1:length(s0_values)
    fprintf('%-10.2f %-10.2f %-10.4f %-10.4f %-10.2f %-10.4f\n', ...
        results.s0(i).s0, results.s0(i).TMRCA, results.s0(i).V_num, results.s0(i).V_an, ...
        results.s0(i).k_final, results.s0(i).f_final);
end

% Сохранение результатов в файл
save('experiment_results.mat', 'results');
fprintf('\nРезультаты сохранены в файл experiment_results.mat\n');

%% ФУНКЦИЯ ДЛЯ ВЫЧИСЛЕНИЯ АНАЛИТИЧЕСКОЙ СКОРОСТИ (V_an)
function V = compute_analytical_velocity(N, s, L, f0, muL)
    % Вычисление аналитической скорости адаптации по формуле
    % V ≈ 2𝑠 log(𝑁√(𝑠𝑈𝑏)) / [log(𝑠/𝑈𝑏 * log(𝑁√(𝑠𝑈𝑏)))]^2
    
    % Частота полезных мутаций на геном
    Ub = muL * (1 - f0);
    
    % Проверка на возможность вычисления
    if Ub <= 0 || s <= 0
        V = 0;
        return;
    end
    
    % Вычисление компонентов формулы
    N_sqrt_sUb = N * sqrt(s * Ub);
    
    % Проверка корректности аргументов логарифмов
    if N_sqrt_sUb <= 1
        V = 0;
        return;
    end
    
    log_N_sqrt_sUb = log(N_sqrt_sUb);
    
    % Второй логарифмический член в знаменателе
    s_Ub_ratio = s / Ub;
    arg_log2 = s_Ub_ratio * log_N_sqrt_sUb;
    
    if arg_log2 <= 1
        V = 0;
        return;
    end
    
    log_arg_log2 = log(arg_log2);
    
    % Вычисление скорости по формуле
    numerator = 2 * s * log_N_sqrt_sUb;
    denominator = log_arg_log2^2;
    
    V = numerator / denominator;
    
    % Дополнительная проверка на физичность
    if V < 0 || V > s * L
        V = 0;
    end
end
