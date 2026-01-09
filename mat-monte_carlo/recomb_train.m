function [TMRCA, adapt] = recomb_train(distribution_s, r, M, s0, L, N, tf, f0, muL, run, exp_name)

global yesgenealogy

%% Аргументы:
% distribution_s - одно из четырех распределений ниже
% r - частота рекомбинаций на геном
% M - число кроссинговеров
% s0 - средний коэффициент отбора
% L - число локусов
% N - численность популяции
% tf - полное время в поколениях
% f0 - начальная частота благоприятного аллеля
% muL - частота мутаций на геном
% run - номер запуска
% exp_name - дополнительное имя эксперимента для названий файлов

%% Включаем подробный вывод для сравнения с C++
debug_output = false;  % Включаем вывод отладочной информации
debug_step = 1;       % Выводим каждый шаг

%% Создание папки для графиков
if ~exist('graphics', 'dir')
    mkdir('graphics');
end

%% Случайный сид для воспроизводимости результатов
rng(run)

%% Частота мутаций на локус
mu = muL/L;                      % временной интервал для одного типа начальных условий
T = 0:tf;                       % времена

%% Распределения s и параметры подписей к рис. 1 

switch distribution_s
    case 'exponential'
        %  g(s) = (1/s0)*exp(-s/s0)
        s = -s0*log(rand(1,L));
        ts = sprintf('экспоненциальное\n exp(-s/s0)\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run);
    case 'const'
        s = s0*ones(1,L);
        ts = sprintf('постоянное s\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run);
    case 'halfgaussian'
        % g(s)=(2/pi/s0)*exp(-s^2/pi/s0^2), s0 avrg
        s = s0*sqrt(pi/2)*abs(randn(1,L));
        ts = sprintf('полугауссово\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run);
end

%% ВЫЧИСЛЕНИЕ АНАЛИТИЧЕСКОЙ СКОРОСТИ АДАПТАЦИИ (V_an)
V_an = compute_analytical_velocity(N, s0, L, f0, muL);

%% Начальные установки

A = (1:N)'*ones(1,L); % матрица предков
tint = round(tf/10); % интервал времени для графиков
col = 'rgbmkrgbmkrgbmkrgbmk';
fsample = 0.1; % процент выборки для пар

%% Начальная популяция

% Случай 1: случайно распределенные благоприятные аллели с фиксированной частотой f0
if f0 ~= 0
    K = (rand(N,L) < f0); 
    % Матрица K: Каждая строка - геном, последовательность 0 (вредные) и 1 (благоприятные)
else
    K = zeros(N,L);
end

% Резервирование памяти
W = zeros(N, length(T));  
P1 = zeros(N, length(T)); 
PL = P1; % Начальные метки родителей для 3 сайтов во времени, для филогении
fsite = zeros(length(T), L);
kav = zeros(1, length(T)); 
Vark = kav; 
fsurvive = kav; 
C = kav; 
Call = kav; 
meanW = kav;
dist_over_L = zeros(1, length(T)); % Добавим для вывода
Knew = zeros(N, L);  % бинарные последовательности ДНК
Anew = zeros(N, L);  % последовательности меток предков
Pnew = zeros(N, L);  % метки родителей

%% ВЫВОД НАЧАЛЬНЫХ УСЛОВИЙ
if debug_output
    fprintf('\n================================================================================\n');
    fprintf('MATLAB DEBUG OUTPUT FOR MONTE CARLO SIMULATION\n');
    fprintf('Parameters: N=%d, L=%d, s0=%.4f, r=%.4f, f0=%.4f, muL=%.4f, tf=%d, M=%d, distribution=%s\n',...
        N, L, s0, r, f0, muL, tf, M, distribution_s);
    fprintf('================================================================================\n');
    
    % Вывод начального состояния
    total_alleles = sum(K(:));
    fprintf('Начальное состояние:\n');
    fprintf('Всего благоприятных аллелей: %d\n', total_alleles);
    fprintf('Среднее на геном: %.4f\n', total_alleles/N);
    fprintf('k_av(0)/L: %.4f\n', total_alleles/N/L);
end

%% Начало эволюции...
for t = T
       
    % симметричные мутации если есть
    if mu > 0
        K = xor(K, rand(N,L) < mu); 
    end
   
    % Начальные метки родителей
    P = (1:N)'*ones(1,L); % матрица родителей
    
    %% Случайная выборка потомков и естественный отбор с "сломанной палкой"
    % столбец N логарифмов приспособленностей; состояние с 0 имеет w=0 (приспособленность 1) по определению
    w = K*(s');
    
    % добавить эпистазис здесь! Входные параметры Tij, Eij устанавливаются отдельно в начале.
    nprogav = exp(w)/mean(exp(w));    % среднее число потомков
    b2 = cumsum(nprogav);
    b1 = [0; b2(1:end-1)];             % сломанная палка
    X = rand(N,1)*N;
    nprog = zeros(N,1);
    for i = 1:N
        nprog(i) = sum(X > b1(i) & X < b2(i)); % актуальное число потомков
    end

    if t <= 10  % Выводим только первые шаги
        fprintf('DEBUG t=%d:\n', t);
        fprintf('  w(1:5) = ');
        fprintf('%.4f ', w(1:5)/s0);
        fprintf('\n');
    
        fprintf('  nprogav = ');
        fprintf('%.4f ', nprogav(1:5));
        fprintf('\n');
        % Выведем nprog для первых 5 особей
        fprintf('  nprog(1:5) = ');
        fprintf('%d ', nprog(1:5));
        fprintf('\n');
    
        % Выведем b1, b2 для первых 5 особей
        fprintf('  b1(1:5) = ');
        fprintf('%.2f ', b1(1:5));
        fprintf('\n');
        fprintf('  b2(1:5) = ');
        fprintf('%.2f ', b2(1:5));
        fprintf('\n');
    
        % Проверим сумму nprog
        fprintf('  sum(nprog) = %d\n', sum(nprog));
    end
  
    %% Обновление популяции
    is = [0; cumsum(nprog(1:(N-1)))];
    for i = 1:N
        if nprog(i)
            Knew(is(i)+1:is(i)+nprog(i),:) = ones(nprog(i),1)*K(i,:); % последовательности ДНК
            Anew(is(i)+1:is(i)+nprog(i),:) = ones(nprog(i),1)*A(i,:); % метки предков
            Pnew(is(i)+1:is(i)+nprog(i),:) = ones(nprog(i),1)*P(i,:); % метки родителей
        end
    end
    K = Knew;
    A = Anew;
    P = Pnew;
     
    %% Рекомбинация случайно выбранных пар с заменой одного родителя
    npairs = round(r*N/2);
    ii = ceil(rand(npairs,2)*N); 
    i1 = ii(:,1); 
    i2 = ii(:,2);  % 2 столбца случайных индексов родителей
    for i = 1:npairs
        % генерация случайных 1 или 0 для каждого сайта, с вероятностями M/L и 1-M/L,
        % соответственно, чтобы отметить кроссинговеры единицами
        % Четные/нечетные xx показывают сегменты сайтов, скопированные от 1-го родителя, 2-го, 1-го, 2-го и т.д.
        xx = cumsum(rand(1,L) < M/L); 
        % сайты, скопированные от 1-го родителя, отмечены как 1
        first = (round(xx/2) == xx/2);     
        % рекомбинантная последовательность ДНК
        prog1 = K(i1(i),:).*first + K(i2(i),:).*(1-first);   
        % рекомбинантная метка предка
        prog1A = A(i1(i),:).*first + A(i2(i),:).*(1-first);
        % рекомбинантная метка родителя
        prog1P = P(i1(i),:).*first + P(i2(i),:).*(1-first);      
        
        %% Замена родителя
        if rand > 0.5  
            K(i1(i),:) = prog1; % ДНК 1-го родителя заменена
            A(i1(i),:) = prog1A; % метка предка 1-го родителя заменена
            P(i1(i),:) = prog1P; % метка предка 1-го родителя заменена
        else
            K(i2(i),:) = prog1; % ДНК 2-го родителя заменена
            A(i2(i),:) = prog1A; % метка предка 2-го родителя заменена
            P(i2(i),:) = prog1P; % метка предка 2-го родителя заменена
        end
    end % спаривания  

    
    %% Запись базовых наблюдаемых величин 
    fsite(t+1,:) = mean(K);        % частоты аллелей в 1-сайте во всех сайтах
    kav(t+1) = L*mean(mean(K));     % среднее число аллелей на геном
    Vark(t+1) = (std(w)/s0)^2;       % дисперсия числа аллелей между геномами
    fsurvive(t+1) = mean(~all(K==0));
    
    % Доля пар гомологичных локусов с общим предком
    xx = round(N*fsample);
    i1 = ceil(N*rand(1,xx)); 
    i2 = ceil(N*rand(1,xx));
    C(t+1) = mean(mean(A(i1,:) == A(i2,:)));
    Call(t+1) = mean(~std(A));
    
    % Приспособленность всех геномов
    W(:,t+1) = w;
    meanW(t+1) = mean(w);
    
    % Вычисление dist (генетическое разнообразие)
    dist_current = 0;
    for locus = 1:L
        f = fsite(t+1, locus);
        dist_current = dist_current + 2.0 * f * (1.0 - f);
    end
    dist_over_L(t+1) = dist_current / L;
    
    %% ПОДРОБНЫЙ ВЫВОД ДЛЯ СРАВНЕНИЯ С C++
    if debug_output && mod(t, debug_step) == 0
        % 2. Статистики по локусам
        min_freq = min(fsite(t+1, :));
        max_freq = max(fsite(t+1, :));
        mean_freq = mean(fsite(t+1, :));
        % 4. Статистики по геномам
        total_alleles = sum(K(:));
        
        % Вывод в табличном формате для удобства
        fprintf('t=%4d: k_av=%.2f, k_av/L=%.4f, V_ark=%.4f, mean_W=%.4f, f_survive=%.4f, C=%.4f, C_all=%.4f, sumK=%d, min_f=%.4f, max_f=%.4f, mean_f=%.4f\n',...
            t, kav(t+1), kav(t+1)/L, Vark(t+1), meanW(t+1), fsurvive(t+1),...
            C(t+1), Call(t+1), total_alleles, min_freq, max_freq, mean_freq);
    end
    
    %% запомнить имена родителей для 1-го и последнего локусов   
    P1(:,t+1) = P(:,1); 
    PL(:,t+1) = P(:,L);    
    
    %% Построение графика волны и спектра ancestral clone в некоторые моменты времени
    if tint*round(t/tint) == t 
        c = col(round(t/tint)+1);
        figure(1)
        subplot(2,2,1)
        [nn,xx] = hist(w);         % гистограмма приспособленности среди геномов
        semilogy(xx/s0, nn, c)
        hold on
        text(xx(end)/s0, nn(end), sprintf('%g',t))
    end % if sampling time

end % цикл по времени
% Эволюция завершена

%% ВЫЧИСЛЕНИЕ ЧИСЛЕННОЙ СКОРОСТИ АДАПТАЦИИ (V_num)
kav = reshape(kav, size(T)); % среднее k

% Скорость адаптации как производная от среднего числа благоприятных аллелей
% Используем центральную разность для внутренних точек
if length(T) > 3
    V_num_vec = zeros(size(kav));
    % Для внутренних точек используем центральную разность
    for i = 2:length(T)-1
        V_num_vec(i) = (kav(i+1) - kav(i-1)) / 2;
    end
    % Для первой и последней точки используем одностороннюю разность
    V_num_vec(1) = kav(2) - kav(1);
    V_num_vec(end) = kav(end) - kav(end-1);
    
    % Средняя численная скорость за последнюю треть времени
    t_start = round(2*tf/3);
    if t_start < length(V_num_vec)
        V_num = mean(V_num_vec(t_start:end));
    else
        V_num = mean(V_num_vec);
    end
else
    % Если мало точек, используем простую разность
    V_num = (kav(end) - kav(1)) / (T(end) - T(1));
    V_num_vec = zeros(size(kav));
    V_num_vec(:) = V_num;
end

%% Финальные графики
figure(1) % Бегущая волна
subplot(2,2,1)
hold off
xlabel('Число аллелей, k');
ylabel('Функция распределения')
title(ts)   % заголовок со значениями параметров
 
%% Средние наблюдаемые величины vs время

Vark = reshape(Vark, size(T)); % Дисперсия k
dist = mean(L*2*fsite.*(1-fsite), 2);
dist = reshape(dist, size(T));

% 1-сайтовая теория
f1site = f0*(f0+(1-f0)*exp(-s0*T)).^(-1); 

%%
subplot(2,2,2)
plot(T, kav/L, 'r', T, sqrt(Vark)/L, 'b', T, dist/L, 'g', T, C, '--k', T, fsurvive, 'k', T, Call, 'm', T, f1site, ':');

% МОДИФИЦИРОВАННЫЙ ЗАГОЛОВОК С ЧИСЛЕННОЙ И АНАЛИТИЧЕСКОЙ СКОРОСТЯМИ
if ~isnan(V_an)
    title(sprintf('V_num=%.3e, V_an=%.3e\n k_{ср}/L кр SD_k/L син, dist/L зел \n C --чер, f_{выж} чер, Call пурп', V_num, V_an));
else
    title(sprintf('V_num=%.3e\n k_{ср}/L кр SD_k/L син, dist/L зел \n C --чер, f_{выж} чер, Call пурп', V_num));
end

xlabel('Время, t');
axis([0 T(end) 0 1])
Cinf = fsurvive(end);

subplot(2,2,3)
meanf_end = mean(fsite(end,:));
for i = 1:L
    plot(T', fsite(:,i), T, f1site, 'k--')
    hold on
end
hold off
ylabel('Частота аллеля во всех локусах')
xlabel('Время, t')
title(sprintf('средн(f_{кон}) = %g', meanf_end))

subplot(2,2,4)
for t = 0:tint:tf
    [nn,xx] = hist(fsite(t+1,:));
    text(xx(end), nn(end), sprintf('t=%g',t))
    plot(xx, nn)
    hold on
end
ylabel('Гистограмма в несколько фиксированных времен')
xlabel('Частота аллеля в локусе, t')
hold off

%% Сохранение Figure 1
filename1 = sprintf('graphics/%s_fig1_run%d_N%d_L%d_s0%.3f_r%.3f', exp_name, run, N, L, s0, r);
saveas(gcf, [filename1, '.png']);

% Вывод скоростей в командное окно
fprintf('Численная скорость адаптации (V_num): %.4f\n', V_num);
fprintf('Аналитическая скорость адаптации (V_an): %.4f\n', V_an);

if yesgenealogy
    %% Траектория генеалогии
    figure(2)
    % Выборка геномов в t=tf
    ii = [1 N/4:N/4:N];
    % Номера генеалогии локусы 1 и 2 
    G1(tf+1,ii) = ii; 
    GL(tf+1,ii) = ii;
    for t = tf:-1:1
        G1(t,ii) = P1(G1(t+1,ii), t+1);
        GL(t,ii) = PL(GL(t+1,ii), t+1);
    end
    subplot(2,2,1)
    plot(T, G1(:,ii)); 
    title(ts)
    ylabel('Родитель, локус 1')

    %% Вычисление TMRCA 
    TMRCA = tf - max(find(~std(G1(:,ii)'))); %#ok<MXFND>
    
    subplot(2,2,3)
    plot(T, GL(:,ii)); 
    xlabel('Время')
    ylabel('Родитель, локус L')

    %% Скорость адаптации за последние 3/4 временного интервала
    t1 = round(tf/4+1); 
    adapt = (meanW(tf+1) - meanW(t1))/(tf-t1);

    %% Траектория приспособленности
    for t = T
        w1(t+1,ii) = W(G1(t+1,ii), t+1);
        wL(t+1,ii) = W(GL(t+1,ii), t+1);
    end
    subplot(2,2,2)
    plot(T, w1)
    xlabel('Время')
    ylabel('Приспособленность')
    subplot(2,2,4)
    plot(T, wL)
    xlabel('Время')
    ylabel('Приспособленность')
    
    %% Сохранение Figure 2
    %filename2 = sprintf('graphics/%s_fig2_run%d_N%d_L%d_s0%.3f_r%.3f', exp_name, run, N, L, s0, r);
    %saveas(gcf, [filename2, '.png']);

    %% Филогенетическое дерево
    figure(3)
    m = length(ii);
    yesconnect1 = ones(size(ii)); % метка построения
    yesconnectL = ones(size(ii)); % метка построения
    color = 'rbmkrbmkrbmk'; % цвета линий

    % Цикл назад во времени 
    for t = tf:-1:0
        % дерево первого сайта
        subplot(2,1,1)
        for i = 1:m 
            % найти все пересечения > i
            jj = find(G1(t+1,ii(i)) == G1(t+1,ii((i+1):m)));
            if isempty(jj) && yesconnect1(i)
                % если нет, построить горизонтальный сегмент
                plot([t t+1], [i i], color(i)); 
                hold on 
            elseif yesconnect1(i)
                % построить наклонный сегмент вверх
                plot([t t+1], [max(jj)+i, i], color(i)); 
                hold on 
                % и больше не строить его
                yesconnect1(i) = 0;  
            end 
        end % по линиям
        
        % дерево последнего сайта
        subplot(2,1,2)
        for i = 1:m 
            % найти все пересечения > i
            jj = find(GL(t+1,ii(i)) == GL(t+1,ii((i+1):m)));
            if isempty(jj) && yesconnectL(i)
                % если нет, построить горизонтальный сегмент
                plot([t t+1], [i i], color(i)); 
                hold on 
            elseif yesconnectL(i)
                % построить наклонный сегмент вверх
                plot([t t+1], [max(jj)+i, i], color(i)); 
                hold on 
                % и больше не строить его
                yesconnectL(i) = 0;  
            end 
        end % по линиям
    end  % по времени

    subplot(2,1,1)
    hold off
    ylabel('Дерево, локус 1')
    xlabel('Время')
    title(ts)
    axis([0 tf+1 0.5 m+0.5])
    subplot(2,1,2)
    hold off
    ylabel('Дерево, локус L')
    xlabel('Время')
    title(ts)
    axis([0 tf+1 0.5 m+0.5])
    
    %% Сохранение Figure 3
    %filename3 = sprintf('graphics/%s_fig3_run%d_N%d_L%d_s0%.3f_r%.3f', exp_name, run, N, L, s0, r);
    %saveas(gcf, [filename3, '.png']);
    
end % yesgenealogy?

fprintf('Графики сохранены в папку graphics с префиксом: %s\n', exp_name);

%% ФУНКЦИИ ДЛЯ ВЫЧИСЛЕНИЯ АНАЛИТИЧЕСКОЙ СКОРОСТИ (V_an)
function V_an = compute_analytical_velocity(N, s, L, f0, muL)
    % Частота полезных мутаций на геном
    Ub = muL * (1 - f0);
    
    % Проверка применимости теории
    if Ub <= 0 || s <= 0
        fprintf('ERROR: Ub or s must be positive\n');
        V_an = 0.0;
        return;
    end
    
    target = log(N);
    
    % Определяем границы поиска V
    V_min = Ub * 1.001;  % V должно быть > Ub для положительного логарифма
    V_max = 1000.0 * s;  % Начальная верхняя граница
    
    % Увеличиваем V_max до тех пор, пока правая часть (51) не станет больше target
    max_iter_expand = 100;
    iter = 0;
    while iter < max_iter_expand && right_side_51(V_max, s, Ub) < target && V_max < 1e20
        V_max = V_max * 2.0;
        iter = iter + 1;
    end
    
    if V_max >= 1e20
        fprintf('ERROR: Cannot find suitable V_max, N might be too large\n');
        V_an = 0.0;
        return;
    end
    
    % Решаем уравнение (51)
    V1 = bisection(V_min, V_max, target, @right_side_51, s, Ub);
    
    % Проверяем условие длинного хвоста
    if V1 * log(V1 / Ub) < s
        fprintf('WARNING: Condition V*ln(V/Ub) >> s might not be satisfied.\n');
        fprintf('  V*ln(V/Ub) = %.4f, s = %.4f\n', V1 * log(V1 / Ub), s);
    end
    
    % Если V1 > s, используем результат (51) для широкого распределения
    if V1 > s
        V_an = V1;
    % Иначе решаем уравнение (52) для узкого распределения
    else
        % Решаем уравнение (52)
        V2 = bisection(V_min, V_max, target, @right_side_52, s, Ub);
        
        % Проверяем условие для формулы (52)
        if V2 < s && V2 > s / log(V2 / Ub)
            V_an = V2;
        else
            % Если условия не выполняются, используем V1 с предупреждением
            fprintf('WARNING: Conditions for formula (52) not satisfied.\n');
            fprintf('  Using V from formula (51): V = %.4f\n', V1);
            V_an = V1;
        end
    end
end

%% Вспомогательные функции
% Правая часть уравнения (51)
function value = right_side_51(V, s, Ub)
    log_term = log(V / (exp(1.0) * Ub));
    term1 = V / (2.0 * s) * (log_term * log_term + 1.0);
    term2 = 0.5 * log((s * s * s * Ub) / (V * V * log(V / Ub)));
    value = term1 - term2;
end

% Правая часть уравнения (52)
function value = right_side_52(V, s, Ub)
    log_term = log(V / (exp(1.0) * Ub));
    term1 = V / (2.0 * s) * (log_term * log_term + 1.0);
    term2 = 0.5 * log((s * s * Ub) / (V * log(V / Ub)));
    value = term1 - term2;
end

% Метод бисекции для решения уравнения f(V) = target
function V_mid = bisection(V_low, V_high, target, func, s, Ub, tol, max_iter)
    % Установка значений по умолчанию для параметров
    if nargin < 8
        max_iter = 1000;
    end
    if nargin < 7
        tol = 1e-8;
    end
    
    f_low = func(V_low, s, Ub) - target;
    f_high = func(V_high, s, Ub) - target;
    
    % Проверка, что функция меняет знак на отрезке
    if f_low * f_high > 0
        % Если оба значения положительны, возможно, V_low слишком велико
        if f_low > 0 && f_high > 0
            V_mid = V_low;
            return;
        end
        % Если оба отрицательны, V_high слишком мало
        V_mid = V_high;
        return;
    end
    
    for iter = 1:max_iter
        V_mid = (V_low + V_high) / 2.0;
        if V_high - V_low < tol
            return;
        end
        
        f_mid = func(V_mid, s, Ub) - target;
        
        if abs(f_mid) < tol
            return;
        elseif f_low * f_mid < 0
            V_high = V_mid;
            f_high = f_mid;
        else
            V_low = V_mid;
            f_low = f_mid;
        end
    end
    
    V_mid = (V_low + V_high) / 2.0;
end
end % конец функции recomb_train