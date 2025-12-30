global yesgenealogy;
yesgenealogy = 1;  % 1 - строить генеалогию, 0 - нет

% Базовые параметры
distribution_s = 'const';
r = 0;      % Частота рекомбинации
M = 3;      % Количество точек рекомбинации
s0_base = 0.1;   % Базовая сила отбора
L_base = 100;    % Длина генома (количество локусов)
N_base = 1000;   % Размер популяции
tf = 150;   % Время моделирования
f0 = 0;     % Начальная частота благоприятных аллелей
muL_base = 0.1;  % Общая частота мутаций
run = 1;    % номер запуска

% Создаем структуру для хранения результатов
results = struct();

%% Серия 4: Влияние силы отбора (s0)
disp('=== СЕРИЯ 4: Влияние силы отбора (s0) ===');
s0_values = [0.05, 0.1, 0.2]; % Меньшая, базовая и большая сила отбора

for i = 1:length(s0_values)
    current_s0 = s0_values(i);
    [TMRCA, adapt_data] = recomb_train(distribution_s, r, M, current_s0, L_base, N_base, tf, f0, muL_base, run, sprintf('exp_s0_%d', i));
    
    % Сохраняем результаты
    results.s0(i).s0 = current_s0;
    results.s0(i).TMRCA = TMRCA;
    results.s0(i).adapt_data = adapt_data;
    
    % Извлекаем данные
    if isstruct(adapt_data) && isfield(adapt_data, 'V_theor')
        results.s0(i).V_theor = adapt_data.V_theor;
    else
        results.s0(i).V_theor = NaN;
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
        results.s0(i).k_av = k_av_data;
        results.s0(i).k_final = k_av_data(end);
        results.s0(i).f_final = k_av_data(end) / L_base;
        
        % Вычисляем фактическую скорость адаптации
        if length(k_av_data) > 1
            time_steps = 1:length(k_av_data);
            results.s0(i).V_actual = mean(diff(k_av_data)) / (time_steps(2)-time_steps(1));
        else
            results.s0(i).V_actual = NaN;
        end
    else
        results.s0(i).k_final = NaN;
        results.s0(i).f_final = NaN;
        results.s0(i).V_actual = NaN;
    end
    
    % Выводим результаты
    fprintf('s0 = %.2f:\n', current_s0);
    fprintf('  TMRCA = %.2f\n', TMRCA);
    if ~isnan(results.s0(i).V_theor)
        fprintf('  Теоретическая скорость адаптации V = %.4f\n', results.s0(i).V_theor);
    end
    if ~isnan(results.s0(i).V_actual)
        fprintf('  Фактическая скорость адаптации V = %.4f\n', results.s0(i).V_actual);
    end
    if ~isnan(results.s0(i).k_final)
        fprintf('  Финальное число благоприятных аллелей k_кон = %.2f\n', results.s0(i).k_final);
        fprintf('  Финальная частота благоприятных аллелей f_кон = %.4f\n\n', results.s0(i).f_final);
    end
end