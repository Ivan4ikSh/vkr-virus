clc;
clear all;

global yesgenealogy;
yesgenealogy = 1;  % 1 - строить генеалогию, 0 - нет

% Количество прогонов для каждого набора параметров
num_runs = 3;

% Базовые параметры (по умолчанию)
base_params.distribution_s = 'const';
base_params.r = 0;      % Частота рекомбинации
base_params.M = 3;      % Количество точек рекомбинации
base_params.s0 = 0.1;   % Базовая сила отбора
base_params.L = 100;    % Длина генома
base_params.N = 1000;   % Размер популяции
base_params.tf = 150;   % Время моделирования
base_params.f0 = 0.01;     % Начальная частота благоприятных аллелей
base_params.muL = 0.01;  % Общая частота мутаций

%% Серия 1: Влияние размера популяции (N)
%{
disp('=== СЕРИЯ 1: Влияние размера популяции (N) ===');
N_values = [1000, 2500, 5000];
exp_params = base_params;

for i = 1:length(N_values)
    exp_params.N = N_values(i);
    results = run_multiple_experiments(exp_params, num_runs, sprintf('s1_N%d', exp_params.N));
    save_results_to_file(results, sprintf('exp1_N_%d_results.txt', exp_params.N));
end

%% Серия 2: Влияние числа локусов (L)
disp('=== СЕРИЯ 2: Влияние числа локусов (L) ===');
L_values = [100, 200, 300];
exp_params = base_params;
exp_params.N = 1000;  % Фиксируем N для серии 2

for i = 1:length(L_values)
    exp_params.L = L_values(i);
    results = run_multiple_experiments(exp_params, num_runs, sprintf('s2_L%d', exp_params.L));
    save_results_to_file(results, sprintf('exp2_L_%d_results.txt', exp_params.L));
end
%}
%% Серия 3: Влияние частоты мутаций (muL)
%{
disp('=== СЕРИЯ 3: Влияние частоты мутаций (muL) ===');
muL_values = [0.01, 0.05];
exp_params = base_params;
exp_params.N = 1000;  % Фиксируем N для серии 3
exp_params.L = 100;   % Фиксируем L для серии 3

for i = 1:length(muL_values)
    exp_params.muL = muL_values(i);
    results = run_multiple_experiments(exp_params, num_runs, sprintf('s3_muL_%.2f', exp_params.muL));
    save_results_to_file(results, sprintf('exp3_muL_%.2f_results.txt', exp_params.muL));
end
%}
%% Серия 4: Влияние силы отбора (s0)
disp('=== СЕРИЯ 4: Влияние силы отбора (s0) ===');
s0_values = [0.05, 0.1];
exp_params = base_params;

for i = 1:length(s0_values)
    exp_params.s0 = s0_values(i);
    results = run_multiple_experiments(exp_params, num_runs, sprintf('s4_s0_%.2f', exp_params.s0));
    save_results_to_file(results, sprintf('exp4_s0_%.2f_results.txt', exp_params.s0));
end

disp('=== ВСЕ ЭКСПЕРИМЕНТЫ ЗАВЕРШЕНЫ ===');
%% Вспомогательные функции

function results = run_multiple_experiments(params, num_runs, exp_tag)
    % Запускает несколько прогонов эксперимента с одними и теми же параметрами
    % и возвращает структуру с усредненными результатами
    
    % Инициализация массивов для хранения результатов всех прогонов
    TMRCA_array = cell(1, num_runs);  % Изменено на cell array для хранения векторов
    V_num_array = zeros(1, num_runs);
    V_an_array = zeros(1, num_runs);
    k_final_array = zeros(1, num_runs);
    f_final_array = zeros(1, num_runs);
    
    fprintf('Эксперимент %s: запуск %d прогонов...\n', exp_tag, num_runs);
    
    for run_num = 1:num_runs
        fprintf('  Прогон %d/%d...\n', run_num, num_runs);
        % Уникальное имя для файла графика
        run_name = sprintf('%s_run%d', exp_tag, run_num);
        % Запуск модели
        [TMRCA, adapt_data] = recomb_train(params.distribution_s, params.r, params.M, params.s0, params.L, params.N, params.tf, params.f0, params.muL, run_num, run_name);
        % Сохраняем TMRCA (может быть вектором)
        TMRCA_array{run_num} = TMRCA;
        % Проверяем структуру adapt_data на наличие нужных полей
        V_num = NaN;
        V_an = NaN;
        
        if isstruct(adapt_data)
            % Извлекаем k_av для расчёта k_кон и f_кон
            if isfield(adapt_data, 'k_av')
                k_av_data = adapt_data.k_av;
                k_final_array(run_num) = k_av_data(end);
                f_final_array(run_num) = k_av_data(end) / params.L;
            else
                % Если нет k_av, ищем числовое поле
                for j = 1:length(field_names)
                    field_data = adapt_data.(field_names{j});
                    if isnumeric(field_data) && length(field_data) > 1
                        k_av_data = field_data;
                        k_final_array(run_num) = k_av_data(end);
                        f_final_array(run_num) = k_av_data(end) / params.L;
                        break;
                    end
                end
            end
        elseif isnumeric(adapt_data)
            % Если adapt_data - это массив (предположительно k_av)
            k_av_data = adapt_data;
            k_final_array(run_num) = k_av_data(end);
            f_final_array(run_num) = k_av_data(end) / params.L;
        else
            k_final_array(run_num) = NaN;
            f_final_array(run_num) = NaN;
        end
        
        % Сохраняем скорости
        V_num_array(run_num) = V_num;
        V_an_array(run_num) = V_an;
    end
    
    % Рассчитываем средние значения и стандартные отклонения
    results = struct();
    
    % Параметры эксперимента
    results.params = params;
    results.exp_tag = exp_tag;
    results.num_runs = num_runs;
    
    % Данные по прогонам
    results.TMRCA_runs = TMRCA_array;  % Сохраняем как cell array
    results.V_num_runs = V_num_array;
    results.V_an_runs = V_an_array;
    results.k_final_runs = k_final_array;
    results.f_final_runs = f_final_array;
    
    % Для TMRCA вычисляем среднее по всем элементам всех векторов
    all_TMRCA_values = [];
    for i = 1:num_runs
        if isnumeric(TMRCA_array{i})
            all_TMRCA_values = [all_TMRCA_values; TMRCA_array{i}(:)];
        end
    end
    
    if ~isempty(all_TMRCA_values)
        results.TMRCA_mean = mean(all_TMRCA_values, 'omitnan');
        results.TMRCA_std = std(all_TMRCA_values, 'omitnan');
    else
        results.TMRCA_mean = NaN;
        results.TMRCA_std = NaN;
    end
    
    % Средние значения и стандартные отклонения для остальных параметров
    results.V_num_mean = mean(V_num_array, 'omitnan');
    results.V_an_mean = mean(V_an_array, 'omitnan');
    results.k_final_mean = mean(k_final_array, 'omitnan');
    results.f_final_mean = mean(f_final_array, 'omitnan');
    
    results.V_num_std = std(V_num_array, 'omitnan');
    results.V_an_std = std(V_an_array, 'omitnan');
    results.k_final_std = std(k_final_array, 'omitnan');
    results.f_final_std = std(f_final_array, 'omitnan');
end

function save_results_to_file(results, filename)
    % Сохраняет результаты эксперимента в текстовый файл
    
    fid = fopen(filename, 'w', 'n', 'UTF-8');
    % Заголовок
    fprintf(fid, 'ЭКСПЕРИМЕНТ: %s\n', results.exp_tag);    
    % Параметры эксперимента
    fprintf(fid, 'ПАРАМЕТРЫ:\n');
    fprintf(fid, 'distribution_s = %s\n', results.params.distribution_s);
    fprintf(fid, 'r = %.2f\n', results.params.r);
    fprintf(fid, 'M = %d\n', results.params.M);
    fprintf(fid, 's0 = %.3f\n', results.params.s0);
    fprintf(fid, 'L = %d\n', results.params.L);
    fprintf(fid, 'N = %d\n', results.params.N);
    fprintf(fid, 'tf = %d\n', results.params.tf);
    fprintf(fid, 'f0 = %.3f\n', results.params.f0);
    fprintf(fid, 'muL = %.3f\n', results.params.muL);
    fprintf(fid, '\n');
    
    % Результаты по прогонам
    fprintf(fid, 'РЕЗУЛЬТАТЫ ПО ПРОГОНАМ:\n');
    fprintf(fid, '----------------------\n');
    fprintf(fid, '%-6s %-15s %-12s %-12s\n', 'Прогон', 'TMRCA', 'k_кон', 'f_кон');
    fprintf(fid, '%-6s %-15s %-12s %-12s\n', '------', '---------------', '------------', '------------');
    
    for i = 1:results.num_runs
        % Форматируем TMRCA для вывода
        if isnumeric(results.TMRCA_runs{i})
            if length(results.TMRCA_runs{i}) == 1
                tmrca_str = sprintf('%.2f', results.TMRCA_runs{i});
            else
                tmrca_str = sprintf('[%dx1]', length(results.TMRCA_runs{i}));
            end
        else
            tmrca_str = 'N/A';
        end
        
        fprintf(fid, '%-6d %-15s %-12.4f %-12.6f\n', i, tmrca_str, results.k_final_runs(i), results.f_final_runs(i));
    end
    fprintf(fid, '\n');
    
    % Статистические показатели
    fprintf(fid, 'СТАТИСТИЧЕСКИЕ ПОКАЗАТЕЛИ (среднее ± стандартное отклонение):\n');
    fprintf(fid, '------------------------------------------------------------\n');
    fprintf(fid, 'TMRCA:                 %.2f ± %.2f поколений\n', results.TMRCA_mean, results.TMRCA_std);
    fprintf(fid, 'Финальное k_кон:       %.4f ± %.4f\n', results.k_final_mean, results.k_final_std);
    fprintf(fid, 'Финальное f_кон:       %.6f ± %.6f\n', results.f_final_mean, results.f_final_std);
    fprintf(fid, '\n');
   
    fprintf('Результаты сохранены в файл: %s\n', filename);
end