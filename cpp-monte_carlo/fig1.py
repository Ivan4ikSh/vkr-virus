import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
import matplotlib.gridspec as gridspec

def load_histogram_data(hist_file):
    """Загрузка данных гистограмм из файла CSV"""
    if not os.path.exists(hist_file):
        print(f"Файл {hist_file} не найден!")
        return None, []
    
    try:
        df = pd.read_csv(hist_file)
        times = df['time'].unique()
        hist_data = {}
        
        for t in times:
            subset = df[df['time'] == t]
            # Сортируем по bin_center для правильного отображения
            subset = subset.sort_values('bin_center')
            hist_data[t] = {
                'bin_centers': subset['bin_center'].values,
                'counts': subset['count'].values
            }
        
        return hist_data, sorted(times)
    except Exception as e:
        print(f"Ошибка загрузки файла {hist_file}: {e}")
        return None, []

def load_subplot2_data(subplot2_file):
    """Загрузка данных для subplot 222"""
    if not os.path.exists(subplot2_file):
        print(f"Файл {subplot2_file} не найден!")
        return None
    
    try:
        df = pd.read_csv(subplot2_file)
        return df
    except Exception as e:
        print(f"Ошибка загрузки файла {subplot2_file}: {e}")
        return None

def load_subplot3_data(subplot3_file):
    """Загрузка данных для subplot 223 (частоты по локусам)"""
    if not os.path.exists(subplot3_file):
        print(f"Файл {subplot3_file} не найден!")
        return None
    
    try:
        df = pd.read_csv(subplot3_file)
        return df
    except Exception as e:
        print(f"Ошибка загрузки файла {subplot3_file}: {e}")
        return None

def load_simulation_params(data_dir):
    """Загрузка параметров симуляции"""
    params_file = os.path.join(data_dir, 'simulation_params.txt')
    if not os.path.exists(params_file):
        print(f"Файл параметров {params_file} не найден!")
        return {}
    
    params = {}
    try:
        with open(params_file, 'r') as f:
            for line in f:
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Преобразуем числовые значения
                    if key in ['N', 'L', 'tf', 'M', 'run']:
                        try:
                            params[key] = int(value)
                        except:
                            params[key] = value
                    elif key in ['s0', 'r', 'f0', 'muL', 'V_an', 'V_num']:
                        try:
                            params[key] = float(value)
                        except:
                            params[key] = value
                    else:
                        params[key] = value
        
        return params
    except Exception as e:
        print(f"Ошибка загрузки параметров: {e}")
        return {}

def plot_figure1_subplot221(hist_data, hist_times, ax, params=None):
    """Построение subplot 221 - гистограммы приспособленности"""
    colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k', 'r', 'g', 'b', 'm', 'c', 'y', 'k']
    
    for i, t in enumerate(hist_times):
        if t not in hist_data:
            continue
            
        data = hist_data[t]
        xx = data['bin_centers']
        nn = data['counts']
        
        # Используем цвет по порядку
        color = colors[i % len(colors)]
        ax.semilogy(xx, nn, color=color, linewidth=1.5, alpha=0.8)
        
        # Добавляем метку времени в конце кривой
        if len(xx) > 0 and len(nn) > 0 and nn[-1] > 0:
            ax.text(xx[-1], nn[-1], f't={t}', fontsize=8, color=color, 
                   ha='left', va='bottom', alpha=0.8)
    
    ax.set_xlabel('Число аллелей, k', fontsize=10)
    ax.set_ylabel('Функция распределения', fontsize=10)
    ax.grid(True, alpha=0.3, which='both', linestyle='--')
    ax.set_ylim(bottom=0.1)
    
    # Формируем заголовок как в MATLAB
    if params:
        dist = params.get('distribution', 'N/A')
        if dist == 'exponential':
            title_text = f'экспоненциальное\n exp(-s/s0)\n'
        elif dist == 'const':
            title_text = f'постоянное s\n'
        elif dist == 'halfgaussian':
            title_text = f'полугауссово\n'
        else:
            title_text = f'{dist}\n'
        
        title_text += (f'N={params.get("N", "N/A")}, r={params.get("r", "N/A")}, '
                      f'L={params.get("L", "N/A")}, s0={params.get("s0", "N/A")}\n'
                      f'f0={params.get("f0", "N/A")}, muL={params.get("muL", "N/A")}, '
                      f'tf={params.get("tf", "N/A")}, M={params.get("M", "N/A")}, '
                      f'run={params.get("run", "N/A")}')
        
        ax.set_title(title_text, fontsize=9, pad=10)

def plot_figure1_subplot222(df, ax, params=None):
    """Построение subplot 222 - основные статистики"""
    if df is None or df.empty:
        ax.text(0.5, 0.5, 'Нет данных для subplot 222', 
                ha='center', va='center', fontsize=12)
        return
    
    # Проверяем наличие необходимых столбцов
    required_columns = ['k_av_L', 'sqrt_Vark_L', 'dist_L', 'C', 'f_survive', 'C_all', 'f1_site']
    available_cols = df.columns.tolist()
    
    # Строим все доступные кривые
    time = df['Time'] if 'Time' in df.columns else range(len(df))
    
    # 1. k_av/L - красный
    if 'k_av_L' in available_cols:
        ax.plot(time, df['k_av_L'], 'r-', linewidth=1.5, label='k_av/L')
    
    # 2. sqrt(Vark)/L - синий
    if 'sqrt_Vark_L' in available_cols:
        ax.plot(time, df['sqrt_Vark_L'], 'b-', linewidth=1.5, label='sqrt(Vark)/L')
    
    # 3. dist/L - зеленый
    if 'dist_L' in available_cols:
        ax.plot(time, df['dist_L'], 'g-', linewidth=1.5, label='dist/L')
    
    # 4. C - черный пунктир
    if 'C' in available_cols:
        ax.plot(time, df['C'], 'k--', linewidth=1.5, label='C')
    
    # 5. f_survive - черный сплошной
    if 'f_survive' in available_cols:
        ax.plot(time, df['f_survive'], 'k-', linewidth=1.5, label='f_survive')
    
    # 6. C_all - пурпурный
    if 'C_all' in available_cols:
        ax.plot(time, df['C_all'], 'm-', linewidth=1.5, label='C_all')
    
    # 7. f1_site - черный с точками
    if 'f1_site' in available_cols:
        ax.plot(time, df['f1_site'], 'k:', linewidth=2.0, label='f1_site')
    
    # Настройка осей
    ax.set_xlabel('Время, t', fontsize=10)
    ax.set_ylabel('Значения', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim([0, time.max() if len(time) > 0 else 1])
    ax.set_ylim([0, 1])
    
    # Заголовок с информацией о скоростях адаптации
    if params:
        v_num = params.get('V_num', 'N/A')
        v_an = params.get('V_an', 'N/A')
        
        if v_num != 'N/A' and v_an != 'N/A':
            title_text = (f'V_num={v_num:.3e}, V_an={v_an:.3e}\n'
                         f'k_ср/L кр SD_k/L син, dist/L зел\n'
                         f'C --чер, f_выж чер, Call пурп')
        else:
            title_text = 'Статистики по времени'
        
        ax.set_title(title_text, fontsize=9, pad=10)
    
    # Легенда (если есть данные)
    if len(ax.lines) > 0:
        ax.legend(loc='best', fontsize=8, framealpha=0.8)

def plot_figure1_subplot223(df, ax, params=None):
    """Построение subplot 223 - частоты аллелей по локусам"""
    if df is None or df.empty:
        ax.text(0.5, 0.5, 'Нет данных для subplot 223', 
                ha='center', va='center', fontsize=12)
        ax.set_xlabel('Время, t', fontsize=10)
        ax.set_ylabel('Частота аллеля во всех локусах', fontsize=10)
        return
    
    # Определяем столбцы с локусами
    locus_columns = [col for col in df.columns 
                     if col.startswith('locus_')]
    
    # Если нет столбцов locus_, ищем альтернативные
    if not locus_columns:
        locus_columns = [col for col in df.columns 
                         if col not in ['Time', 'f1_site_theoretical'] 
                         and not col.startswith('Unnamed')]
    
    if not locus_columns:
        ax.text(0.5, 0.5, 'Нет данных о локусах', 
                ha='center', va='center', fontsize=12)
        return
    
    time = df['Time'] if 'Time' in df.columns else range(len(df))
    
    # Рисуем частоты для каждого локуса (тонкие линии)
    for locus_col in locus_columns:
        ax.plot(time, df[locus_col], linewidth=0.5, alpha=0.6)
    
    # Рисуем теоретическую частоту f1_site (черный пунктир)
    if 'f1_site_theoretical' in df.columns:
        ax.plot(time, df['f1_site_theoretical'], 'k--', 
                linewidth=2.0, label='Теоретическая f1_site', alpha=0.8)
    
    # Вычисляем среднюю частоту в последний момент времени
    if len(df) > 0 and len(locus_columns) > 0:
        last_row = df.iloc[-1]
        mean_f_end = last_row[locus_columns].mean()
        title_suffix = f'средн(f_кон) = {mean_f_end:.4f}'
    else:
        title_suffix = ''
    
    ax.set_xlabel('Время, t', fontsize=10)
    ax.set_ylabel('Частота аллеля во всех локусах', fontsize=10)
    ax.set_title(f'Частоты аллелей во всех локусах\n{title_suffix}', 
                 fontsize=9, pad=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim([0, time.max() if len(time) > 0 else 1])
    ax.set_ylim([0, 1])
    
    # Легенда только для теоретической кривой
    if 'f1_site_theoretical' in df.columns:
        ax.legend(loc='best', fontsize=8, framealpha=0.8)

def plot_figure1_subplot224(hist_data, hist_times, ax, params=None):
    """Построение subplot 224 - гистограммы частот аллелей"""
    if not hist_data:
        ax.text(0.5, 0.5, 'Нет данных для subplot 224', 
                ha='center', va='center', fontsize=12)
        return
    
    colors = ['r', 'g', 'b', 'm', 'c', 'y', 'k', 'r', 'g', 'b', 'm', 'c', 'y', 'k']
    
    for i, t in enumerate(hist_times):
        if t not in hist_data:
            continue
            
        data = hist_data[t]
        xx = data['bin_centers']
        nn = data['counts']
        
        # Используем цвет по порядку (как в MATLAB)
        color = colors[i % len(colors)]
        
        # Рисуем линейный график (как в MATLAB, не semilogy!)
        ax.plot(xx, nn, color=color, linewidth=1.5, alpha=0.8)
        
        # Добавляем метку времени в конце кривой
        if len(xx) > 0 and len(nn) > 0:
            # Находим последнюю точку с ненулевым значением
            valid_indices = np.where(nn > 0)[0]
            if len(valid_indices) > 0:
                last_idx = valid_indices[-1]
                ax.text(xx[last_idx], nn[last_idx], f't={t}', 
                       fontsize=8, color=color, ha='left', va='bottom', alpha=0.8)
    
    ax.set_xlabel('Частота аллеля в локусе', fontsize=10)
    ax.set_ylabel('Гистограмма', fontsize=10)
    ax.set_title('Гистограмма в несколько фиксированных времен', 
                 fontsize=9, pad=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Устанавливаем разумные пределы по осям
    if hist_data:
        all_xx = []
        all_nn = []
        for t in hist_times:
            if t in hist_data:
                all_xx.extend(hist_data[t]['bin_centers'])
                all_nn.extend(hist_data[t]['counts'])
        
        if all_xx:
            ax.set_xlim([min(all_xx) * 0.9, max(all_xx) * 1.1])
        if all_nn:
            ax.set_ylim([0, max(all_nn) * 1.1])

def plot_figure1_all(data_dir, output_file=None):
    """Построение всех 4 subplots (как в MATLAB)"""
    # Загружаем параметры симуляции
    params = load_simulation_params(data_dir)
    
    # Загружаем данные для каждого subplot
    # Subplot 221: гистограммы приспособленности
    hist_file_221 = os.path.join(data_dir, 'fig1_221.csv')
    hist_data_221, hist_times_221 = load_histogram_data(hist_file_221)
    
    # Subplot 222: основные статистики
    subplot2_file = os.path.join(data_dir, 'fig1_222.csv')
    df_subplot2 = load_subplot2_data(subplot2_file)
    
    # Subplot 223: частоты по локусам
    subplot3_file = os.path.join(data_dir, 'fig1_223.csv')
    df_subplot3 = load_subplot3_data(subplot3_file)
    
    # Subplot 224: гистограммы частот аллелей
    hist_file_224 = os.path.join(data_dir, 'fig1_224.csv')
    hist_data_224, hist_times_224 = load_histogram_data(hist_file_224)
    
    # Создаем фигуру с 4 subplots (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax1, ax2, ax3, ax4 = axes.flatten()
    
    # Subplot 221: гистограммы приспособленности
    if hist_data_221 is not None and hist_times_221:
        plot_figure1_subplot221(hist_data_221, hist_times_221, ax1, params)
    else:
        ax1.text(0.5, 0.5, 'Нет данных для subplot 221', 
                ha='center', va='center', fontsize=14)
    
    # Subplot 222: основные статистики
    plot_figure1_subplot222(df_subplot2, ax2, params)
    
    # Subplot 223: частоты по локусам
    plot_figure1_subplot223(df_subplot3, ax3, params)
    
    # Subplot 224: гистограммы частот аллелей
    plot_figure1_subplot224(hist_data_224, hist_times_224, ax4, params)
    
    # Общий заголовок для всей фигуры (опционально)
    if params:
        exp_name = params.get('exp_name', 'experiment')
        fig.suptitle(f'Результаты симуляции: {exp_name}', fontsize=12, y=0.98)
    
    plt.tight_layout()
    
    # Сохраняем график
    if output_file:
        # Создаем директорию, если нужно
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Полная figure1 сохранена в: {output_file}")
    
    plt.show()
    
    # Выводим информацию о загруженных данных
    print("\nСтатус загрузки данных:")
    print(f"  Subplot 221: {'✓' if hist_data_221 else '✗'} {hist_file_221}")
    print(f"  Subplot 222: {'✓' if df_subplot2 is not None else '✗'} {subplot2_file}")
    print(f"  Subplot 223: {'✓' if df_subplot3 is not None else '✗'} {subplot3_file}")
    print(f"  Subplot 224: {'✓' if hist_data_224 else '✗'} {hist_file_224}")

def main():
    parser = argparse.ArgumentParser(
        description='Построение figure1 для симуляции популяционной генетики',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Пример использования:
  python fig1.py --dir results/experiment1 --output plots/figure1.png""")
    
    parser.add_argument('--dir', type=str, required=True, 
                       help='Директория с результатами симуляции')
    parser.add_argument('--output', type=str, default=None, 
                       help='Имя файла для сохранения графика (например, figure1.png)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.dir):
        print(f"Ошибка: директория {args.dir} не существует!")
        return
    
    print(f"Загружаем данные из директории: {args.dir}")
    plot_figure1_all(args.dir, args.output)
    print("Готово!")

if __name__ == "__main__":
    main()