import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import argparse
import os
from pathlib import Path

# Настройка для поддержки русского языка
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
matplotlib.rcParams['font.family'] = 'sans-serif'

def load_parameters(parameters_file):
    """Загрузка параметров из файла"""
    params = {}
    with open(parameters_file, 'r') as f:
        for line in f:
            if '=' in line:
                key, value = line.strip().split('=')
                try:
                    params[key] = float(value) if '.' in value else int(value)
                except:
                    params[key] = value
    return params

def plot_figure1(data_dir):
    """Построение первого графика (аналог MATLAB fig1)"""
    
    # Загружаем параметры
    params_file = os.path.join(data_dir, 'parameters.txt')
    params = load_parameters(params_file)
    
    # Создаем фигуру с 4 подграфиками
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Подграфик 1: Гистограммы распределения приспособленностей
    print("Загружаем гистограммы приспособленностей...")
    fitness_hist_file = os.path.join(data_dir, 'fitness_histograms.csv')
    if os.path.exists(fitness_hist_file):
        fitness_data = pd.read_csv(fitness_hist_file)
        
        # Получаем уникальные времена
        times = fitness_data['time'].unique()
        
        # Цвета для разных времен
        colors = plt.cm.tab10(np.linspace(0, 1, len(times)))
        
        for t, color in zip(times, colors):
            t_data = fitness_data[fitness_data['time'] == t]
            axes[0, 0].semilogy(t_data['bin_center'] / params['s0'], t_data['count'], color=color, linewidth=1)
            # Добавляем метку времени
            if len(t_data) > 0:
                axes[0, 0].text(t_data['bin_center'].iloc[-1] / params['s0'], t_data['count'].iloc[-1], f't={t}', fontsize=8, color=color)

        axes[0, 0].set_xlabel('Число аллелей, k (нормировано на s0)')
        axes[0, 0].set_ylabel('Функция распределения (log)')
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].set_title('Распределение приспособленностей в разные моменты времени')
    
    # 2. Подграфик 2: Средние наблюдаемые величины
    print("Загружаем средние наблюдаемые величины...")
    mean_file = os.path.join(data_dir, 'mean_observables.csv')
    if os.path.exists(mean_file):
        mean_data = pd.read_csv(mean_file)
        T = mean_data['t'].values
        
        axes[0, 1].plot(T, mean_data['kav_over_L'], 'r-', linewidth=2, label=r'$k_{ср}/L$')
        axes[0, 1].plot(T, mean_data['sqrtVark_over_L'], 'b-', linewidth=2, label=r'$\sqrt{Var_k}/L$')
        axes[0, 1].plot(T, mean_data['dist_over_L'], 'g-', linewidth=2, label=r'$dist/L$')
        axes[0, 1].plot(T, mean_data['C'], 'k--', linewidth=1.5, label='C')
        axes[0, 1].plot(T, mean_data['fsurvive'], 'k-', linewidth=1.5, label=r'$f_{выж}$')
        axes[0, 1].plot(T, mean_data['Call'], 'm-', linewidth=1.5, label='Call')
        axes[0, 1].plot(T, mean_data['f1site'], ':', linewidth=2, label=r'$f_{1сайт}$ теория')
        
        # Добавляем заголовок со скоростями
        title_text = (f'Сравнение численной и аналитической скоростей адаптации\n'
                     f'V_num={params.get("V_num", 0):.3e}, V_an={params.get("V_an", 0):.3e}\n'
                     f'k_ср/L (кр), SD_k/L (син), dist/L (зел)\n'
                     f'C (--чер), f_выж (чер), Call (пурп)')
        axes[0, 1].set_title(title_text)
        axes[0, 1].set_xlabel('Время, t')
        axes[0, 1].set_ylabel('Нормированные величины')
        axes[0, 1].legend(loc='best', fontsize=9)
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_ylim([0, 1])
    
    # 3. Подграфик 3: Частоты аллелей во всех локусах
    print("Загружаем частоты аллелей по локусам...")
    fsite_file = os.path.join(data_dir, 'fsite_matrix.csv')
    if os.path.exists(fsite_file):
        # Загружаем все данные
        fsite_data = pd.read_csv(fsite_file)
        
        # Для каждого локуса строим график
        loci = fsite_data['locus'].unique()
        # Берем только первые 20 локусов для наглядности
        if len(loci) > 20:
            loci = loci[:20]
        
        for locus in loci:
            locus_data = fsite_data[fsite_data['locus'] == locus]
            axes[1, 0].plot(locus_data['time'], locus_data['frequency'], 
                          alpha=0.6, linewidth=0.8)
        
        # Добавляем теоретическую кривую
        if 'f1site' in mean_data.columns:
            axes[1, 0].plot(T, mean_data['f1site'], 'k--', linewidth=2, 
                          label='Теория (1 сайт)')
        
        # Вычисляем среднюю конечную частоту
        last_time = fsite_data['time'].max()
        last_freqs = fsite_data[fsite_data['time'] == last_time]['frequency']
        mean_f_end = last_freqs.mean()
        
        axes[1, 0].set_title(f'Частота аллеля во всех локусах\nсредн(f_кон) = {mean_f_end:.3f}')
        axes[1, 0].set_xlabel('Время, t')
        axes[1, 0].set_ylabel('Частота аллеля')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].legend()
    
    # 4. Подграфик 4: Гистограмма частот аллелей в разные моменты времени
    print("Загружаем гистограммы частот...")
    freq_hist_file = os.path.join(data_dir, 'frequency_histograms.csv')
    if os.path.exists(freq_hist_file):
        freq_hist_data = pd.read_csv(freq_hist_file)
        
        # Получаем уникальные времена
        freq_times = freq_hist_data['time'].unique()
        
        for t in freq_times:
            t_data = freq_hist_data[freq_hist_data['time'] == t]
            axes[1, 1].plot(t_data['bin_center'], t_data['count'], linewidth=1.5, label=f't={t}')
            # Добавляем метку времени
            if len(t_data) > 0:
                axes[1, 1].text(t_data['bin_center'].iloc[-1], t_data['count'].iloc[-1], f't={t}', fontsize=8)
        
        axes[1, 1].set_title('Гистограмма частот аллелей в разные моменты времени')
        axes[1, 1].set_xlabel('Частота аллеля в локусе')
        axes[1, 1].set_ylabel('Количество локусов')
        axes[1, 1].grid(True, alpha=0.3)
    
    # Добавляем общий заголовок с параметрами
    dist_str = params.get('distribution', 'const')
    if dist_str == 'const':
        dist_name = 'постоянное'
    elif dist_str == 'exponential':
        dist_name = 'экспоненциальное'
    else:
        dist_name = 'полугауссово'
    
    suptitle = (f'Распределение s: {dist_name}\n'
               f'N={params.get("N", 0)}, r={params.get("r", 0)}, L={params.get("L", 0)}, '
               f's0={params.get("s0", 0)}\n'
               f'f0={params.get("f0", 0)}, muL={params.get("muL", 0)}, '
               f'tf={params.get("tf", 0)}, M={params.get("M", 0)}, run={params.get("run", 0)}')
    
    plt.suptitle(suptitle, fontsize=12, y=1.02)
    
    # Настраиваем layout
    plt.tight_layout()
    
    # Сохраняем график
    output_file = os.path.join(data_dir, 'figure1.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"График сохранен: {output_file}")
    
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Построение графиков для симуляции популяционной генетики')
    parser.add_argument('--dir', type=str, required=True, help='Директория с данными для построения графиков')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.dir):
        print(f"Ошибка: директория {args.dir} не существует!")
        return
    
    print(f"Построение графиков для данных из: {args.dir}")
    plot_figure1(args.dir)

if __name__ == "__main__":
    main()