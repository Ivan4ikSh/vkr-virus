% --- Высокопроизводительная графовая модель антигенной эволюции ---
clear; clc;
%% Model parameters
R0 = 1.8;                  % basic reproduction ratio
Ub = 1e-3;                 % mutation rate
a = 7;                     % immunity half-distance
N = 1e8;                   % total population size
%% Internal parameters
dt = 0.5;                  % time step
M = 3000;      % number of time points
SEED = 1;
% Фиксируем сид для воспроизводимости результатов
rng(SEED); 

% 2. Выделение памяти
max_nodes = 8000; % Запас памяти под узлы
I = zeros(1, max_nodes); 
R = zeros(1, max_nodes); 
I(1) = 1 / N; % Стартовый посев: 100 больных человек
R(1) = 0.0;
D = zeros(max_nodes, max_nodes);
K = zeros(max_nodes, max_nodes);
K(1,1) = 1; 
sources = []; targets = [];
node_counter = 1;
active_nodes = 1; 
primary_child = zeros(1, max_nodes); % Память об основном потомке штамма
history_I = zeros(M, max_nodes);
fprintf('=== Запуск симуляции на %d итераций ===\n', M);

% 3. Главный цикл
for t = 1:M
    % --- Отладочный вывод с нормами ---
    if mod(t, round(M/10)) == 0
        sum_I = sum(I(active_nodes));
        sum_R = sum(R(1:node_counter));
        fprintf('Шаг %d/%d | Узлов: %d | Активных: %d | Сумма I: %.4f | Сумма R: %.2f\n', ...
                t, M, node_counter, length(active_nodes), sum_I, sum_R);
    end
    
    % --- А. Отсечение мертвых ветвей (Pruning) ---
    new_active = [];
    % Стохастическая функция сама обнуляет I, поэтому порог строго > 0
    for idx = 1:length(active_nodes)
        i = active_nodes(idx);
        if I(i) > 0
            new_active(end+1) = i; 
        end
    end
    active_nodes = new_active; 
    
    if isempty(active_nodes)
        fprintf('ВНИМАНИЕ: Инфекция вымерла на шаге %d!\n', t);
        break; 
    end
    
    new_I = I;
    new_R = R;
    nodes_to_add = [];
    
    % --- Б. Расчет SIR ---
    K_active_Q = K(active_nodes, 1:node_counter);
    R_all = R(1:node_counter)';
    
    K_active_P = K(active_nodes, active_nodes);
    I_active = I(active_nodes)';
    
    Q_vec = 1 - (K_active_Q * R_all); 
    P_vec = K_active_P * I_active;
    
    for idx = 1:length(active_nodes)
        i = active_nodes(idx);
        
        Q_i = Q_vec(idx);
        P_i = P_vec(idx);
        
        % 1. Детерминированный шаг SIR
        new_I(i) = I(i) * (1 + dt * (R0 * Q_i - 1));
        new_R(i) = R(i) * (1 - dt * R0 * P_i) + dt * I(i);
        
        new_I(i) = max(0, new_I(i));
        new_R(i) = max(0, new_R(i));
        
        % 2. СТОХАСТИЧЕСКИЙ ДРЕЙФ для текущего штамма (Broken stick)
        new_I(i) = apply_broken_stick(new_I(i), N);
        
        % Если штамм вымер после дрейфа, нет смысла считать для него мутации
        if new_I(i) == 0
            continue;
        end
    end
        
        % 3. МУТАЦИИ с применением Broken stick
    expected_mutant_density = Ub * I(i) * dt; 
    actual_mutant_density = apply_broken_stick(expected_mutant_density, N);

    if actual_mutant_density > 0
        % Логика группировки мутаций (спасение от комбинаторного взрыва)
        % С вероятностью 95% мутанты сливаются в основную ветвь (как в 2D сетке).
        % С вероятностью 5% происходит "видообразование" - новая параллельная ветвь.
        is_divergence = rand() < 0.05; 

        if primary_child(i) == 0 || is_divergence
            % --- А. Создаем НОВУЮ ветвь ---
            node_counter = node_counter + 1;
            child_node = node_counter;

            % Запоминаем первого потомка как "основного"
            if primary_child(i) == 0
                primary_child(i) = child_node;
            end

            if child_node > size(D, 1)
                D(child_node+1000, child_node+1000) = 0;
                K(child_node+1000, child_node+1000) = 0;
            end

            sources(end+1) = i;
            targets(end+1) = child_node;

            D(child_node, 1:node_counter-1) = D(i, 1:node_counter-1) + 1;
            D(1:node_counter-1, child_node) = D(1:node_counter-1, i) + 1;
            D(child_node, child_node) = 0;

            K(child_node, 1:node_counter-1) = exp(-D(child_node, 1:node_counter-1) / a);
            K(1:node_counter-1, child_node) = exp(-D(1:node_counter-1, child_node) / a);
            K(child_node, child_node) = 1;

            new_I(child_node) = actual_mutant_density;
            new_I(i) = max(0, new_I(i) - actual_mutant_density);
            new_R(child_node) = 0;

            nodes_to_add(end+1) = child_node;

        else
            % --- Б. Слияние с СУЩЕСТВУЮЩИМ потомком ---
            child_node = primary_child(i);
            new_I(child_node) = new_I(child_node) + actual_mutant_density;
            new_I(i) = max(0, new_I(i) - actual_mutant_density);

            % Если потомок спал, "будим" его
            if new_I(child_node) > 0 && ~ismember(child_node, active_nodes) && ~ismember(child_node, nodes_to_add)
                nodes_to_add(end+1) = child_node;
            end
        end
    end

    I(1:node_counter) = new_I(1:node_counter);
    R(1:node_counter) = new_R(1:node_counter);
    active_nodes = [active_nodes, nodes_to_add];
    
    history_I(t, 1:node_counter) = I(1:node_counter);
end
fprintf('Симуляция завершена. Сохранение графиков...\n');

% Подготовка данных
final_nodes = min(node_counter, size(D, 1));
history_I = history_I(1:t, 1:final_nodes);
D_final = D(1:final_nodes, 1:final_nodes);
K_final = K(1:final_nodes, 1:final_nodes);

% =========================================================
% ВИЗУАЛИЗАЦИЯ И ЭКСПОРТ (4 окна)
% =========================================================

% Формируем строку с параметрами для добавления в заголовки
param_str = sprintf('(R0=%.1f, Ub=%g, a=%.1f)', R0, Ub, a);

% --- Окно 1: Динамика волн ---
fig1 = figure('Name', 'Плотность инфицированных (I)', 'Color', 'w', 'Position', [50, 400, 700, 350]);
plot_data = history_I;
plot_data(plot_data == 0) = NaN; % Заменяем нули на NaN для красивой лог-шкалы

semilogy(plot_data, 'LineWidth', 1.0);
title(sprintf('Эволюция вирусных линий %s\nШагов: %d, Узлов: %d', param_str, t, final_nodes));
xlabel('Время (итерации)'); 
ylabel('Доля инфицированных (log)'); 
grid on;
ylim([0.5 / N, max(history_I(:)) * 2]);
saveas(fig1, 'Fig1_Infection_Dynamics.png');

% --- Окно 2: Антигенное дерево ---
fig2 = figure('Name', 'Антигенное дерево', 'Color', 'w', 'Position', [800, 400, 600, 600]);

% Строим граф с жесткой привязкой индексов (чтобы гарантировать порядок узлов)
final_tree = digraph();
final_tree = addnode(final_tree, final_nodes);
if ~isempty(sources)
    valid_edges = (sources <= final_nodes) & (targets <= final_nodes);
    final_tree = addedge(final_tree, sources(valid_edges), targets(valid_edges));
end

% Отрисовываем и СОХРАНЯЕМ КООРДИНАТЫ для 4-го окна!
p_main = plot(final_tree, 'Layout', 'layered', 'Direction', 'down');
tree_x = p_main.XData; 
tree_y = p_main.YData;

% Ищем активные и вымершие штаммы на финальный момент
active_at_end = find(history_I(t, 1:final_nodes) > 0);
extinct_at_end = find(history_I(t, 1:final_nodes) == 0);

p_main.EdgeColor = [0.8 0.8 0.8];

% Зеленые КРУГИ - история (вымершие)
highlight(p_main, extinct_at_end, 'NodeColor', [0.4 0.8 0.4], 'Marker', 'o', 'MarkerSize', 5);
% КРАСНЫЕ КРУГИ - передовой край (активные)
highlight(p_main, active_at_end, 'NodeColor', [1 0 0], 'Marker', 'o', 'MarkerSize', 8);

title(sprintf('Филогенетическое дерево %s\n(Красные - активные, Зеленые - вымершие)', param_str)); 
set(gca, 'Visible', 'off');
saveas(fig2, 'Fig2_Phylogenetic_Tree.png');

% --- Окно 3: Матрицы D и K ---
fig3 = figure('Name', 'Анализ антигенного пространства', 'Color', 'w', 'Position', [150, 50, 900, 400]);
subplot(1, 2, 1);
imagesc(D_final); colorbar; colormap(gca, 'parula'); axis square;
title(sprintf('Матрица дистанций D\n%s', param_str)); xlabel('Индекс штамма'); ylabel('Индекс штамма');

subplot(1, 2, 2);
imagesc(K_final); colorbar; colormap(gca, 'hot'); axis square;
title(sprintf('Матрица кросс-иммунитета K\n%s', param_str)); xlabel('Индекс штамма'); ylabel('Индекс штамма');
saveas(fig3, 'Fig3_Antigenic_Matrices.png');

% --- Окно 4: 8 подграфиков (Срезы плотности + Рост дерева) ---
fig4 = figure('Name', 'Поперечные срезы и рост дерева', 'Color', 'w', 'Position', [50, 50, 1400, 600]);

% Вычисляем 4 точки во времени
time_points = round([M/4, M/2, 3*M/4, M]);

for idx_t = 1:4
    t_snap = min(time_points(idx_t), t); 
    
    % --- ВЕРХНИЙ РЯД: Срез плотности I (subplot 1-4) ---
    subplot(2, 4, idx_t);
    
    I_snap = history_I(t_snap, 1:final_nodes);
    I_snap(I_snap == 0) = NaN; 
    
    semilogy(I_snap, 'LineWidth', 1.5, 'Color', 'b');
    
    title(sprintf('Шаг %d\n%s', t_snap, param_str), 'FontSize', 9);
    xlabel('Индекс штамма');
    if idx_t == 1
        ylabel('Доля инфицированных (log)'); 
    end
    grid on;
    ylim([0.5 / N, max(history_I(:)) * 2]); 
    xlim([1, final_nodes]); 
    
    % --- НИЖНИЙ РЯД: Состояние дерева в этот же момент (subplot 5-8) ---
    subplot(2, 4, idx_t + 4);
    
    % Определяем, какие узлы существовали на момент t_snap
    active_anytime = any(history_I(1:t_snap, 1:final_nodes) > 0, 1);
    N_snap = find(active_anytime, 1, 'last');
    
    if isempty(N_snap)
        N_snap = 1; 
    end
    
    % Строим граф только до N_snap
    snap_tree = digraph();
    snap_tree = addnode(snap_tree, N_snap);
    valid_edges_snap = (sources <= N_snap) & (targets <= N_snap);
    if any(valid_edges_snap)
        snap_tree = addedge(snap_tree, sources(valid_edges_snap), targets(valid_edges_snap));
    end
    
    % Рисуем дерево, ЖЕСТКО задавая координаты из финального дерева
    p_snap = plot(snap_tree, 'XData', tree_x(1:N_snap), 'YData', tree_y(1:N_snap));
    p_snap.EdgeColor = [0.8 0.8 0.8];
    
    % Красим узлы на момент t_snap
    active_snap = find(history_I(t_snap, 1:N_snap) > 0);
    extinct_snap = find(history_I(t_snap, 1:N_snap) == 0);
    
    highlight(p_snap, extinct_snap, 'NodeColor', [0.4 0.8 0.4], 'Marker', 'o', 'MarkerSize', 4);
    highlight(p_snap, active_snap, 'NodeColor', [1 0 0], 'Marker', 'o', 'MarkerSize', 7);
    
    title(sprintf('Рост дерева (до узла %d)', N_snap), 'FontSize', 9);
    set(gca, 'Visible', 'off');
    
    % Фиксируем границы осей, чтобы рамка не скакала
    xlim([min(tree_x)-0.5, max(tree_x)+0.5]);
    ylim([min(tree_y)-0.5, max(tree_y)+0.5]);
end
saveas(fig4, 'Fig4_Snapshots_and_Trees.png');

fprintf('Все графики успешно сохранены в текущую папку как PNG файлы.\n');
% =========================================================
% ЛОКАЛЬНЫЕ ФУНКЦИИ
% =========================================================
function out = apply_broken_stick(in, N)
    % Функция реализует стохастический демографический дрейф
    % (аппроксимация распределения Пуассона для малых чисел)
    exp_num = in * N;
    if exp_num <= 10 && exp_num > 0
        xm = round(6 * exp_num);
        if xm > 0
            out = sum(xm * rand(1, xm) < exp_num) / N;
        else
            out = 0;
        end
    elseif exp_num <= 0
        out = 0;
    else
        % Для больших чисел ( > 10 ) оставляем детерминированную динамику
        out = in;
    end
end