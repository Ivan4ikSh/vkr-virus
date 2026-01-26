import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.colors as mcolors
import os
import glob
import re
import pandas as pd
from matplotlib import cm

class VirusModelVisualizer:
    def __init__(self, output_dir="output"):
        """
        Инициализация визуализатора для модели вирусной эволюции
        
        Parameters:
        -----------
        output_dir : str
            Папка с результатами симуляции
        """
        self.output_dir = output_dir
        self.fig = None
        self.axes = None
        self.anim = None
        
        # Проверяем существование папки
        if not os.path.exists(output_dir):
            print(f"Ошибка: Папка {output_dir} не существует!")
            return
            
        # Находим все файлы с результатами
        self.find_data_files()
        
    def find_data_files(self):
        """Находит все файлы с результатами в папке output"""
        # Файлы с матрицами на каждом шаге
        self.I_files = sorted(glob.glob(os.path.join(self.output_dir, "state_I_step_*.csv")), key=self.extract_step_number)
        self.S_files = sorted(glob.glob(os.path.join(self.output_dir, "state_S_step_*.csv")), key=self.extract_step_number)
        self.R_files = sorted(glob.glob(os.path.join(self.output_dir, "state_R_step_*.csv")), key=self.extract_step_number)
        
        # Финальные файлы
        self.final_I_file = os.path.join(self.output_dir, "final_I.csv")
        self.final_S_file = os.path.join(self.output_dir, "final_S.csv")
        self.final_R_file = os.path.join(self.output_dir, "final_R.csv")
        
        # Временные ряды
        self.norm_file = os.path.join(self.output_dir, "norm_time_series.csv")
        self.finf_file = os.path.join(self.output_dir, "finf_time_series.csv")
        
        print(f"Найдено {len(self.I_files)} шагов симуляции")
        
    def extract_step_number(self, filename):
        """Извлекает номер шага из имени файла"""
        match = re.search(r'step_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    def load_matrix(self, filename):
        """Загружает матрицу из CSV файла"""
        return np.loadtxt(filename, delimiter=',')
    
    def load_time_series(self):
        """Загружает временные ряды"""
        norm_data = pd.read_csv(self.norm_file) if os.path.exists(self.norm_file) else None
        finf_data = pd.read_csv(self.finf_file) if os.path.exists(self.finf_file) else None
        return norm_data, finf_data
    
    def create_composite_plot(self, step_idx=0):
        """
        Создает составной график для конкретного шага
        
        Parameters:
        -----------
        step_idx : int
            Индекс шага (по порядку в списке файлов)
        """
        if step_idx >= len(self.I_files):
            print(f"Ошибка: шаг {step_idx} не существует!")
            return
            
        # Загружаем данные для данного шага
        I = self.load_matrix(self.I_files[step_idx])
        S = self.load_matrix(self.S_files[step_idx])
        R = self.load_matrix(self.R_files[step_idx])
        
        step_num = self.extract_step_number(self.I_files[step_idx])
        
        # Создаем фигуру с 4 subplots
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Настраиваем общий заголовок
        fig.suptitle(f'Двумерная модель антигенной эволюции вируса\nШаг: {step_num}', fontsize=16)
        
        # 1. Инфицированные (I) - heatmap
        im1 = axes[0, 0].imshow(I, cmap='hot', aspect='auto', origin='lower')
        axes[0, 0].set_title('Инфицированные (I)')
        axes[0, 0].set_xlabel('Антигенная координата (x)')
        axes[0, 0].set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im1, ax=axes[0, 0], label='Доля инфицированных')
        
        # 2. Восприимчивые (S) - heatmap
        im2 = axes[0, 1].imshow(S, cmap='Blues', aspect='auto', origin='lower')
        axes[0, 1].set_title('Восприимчивые (S)')
        axes[0, 1].set_xlabel('Антигенная координата (x)')
        axes[0, 1].set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im2, ax=axes[0, 1], label='Доля восприимчивых')
        
        # 3. Выздоровевшие (R) - heatmap
        im3 = axes[1, 0].imshow(R, cmap='Greens', aspect='auto', origin='lower')
        axes[1, 0].set_title('Выздоровевшие (R)')
        axes[1, 0].set_xlabel('Антигенная координата (x)')
        axes[1, 1].set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im3, ax=axes[1, 0], label='Доля выздоровевших')
        
        # 4. Суммарное распределение (I+S+R) - heatmap
        total = I + S + R
        im4 = axes[1, 1].imshow(total, cmap='viridis', aspect='auto', origin='lower')
        axes[1, 1].set_title('Суммарное распределение (I+S+R)')
        axes[1, 1].set_xlabel('Антигенная координата (x)')
        axes[1, 1].set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im4, ax=axes[1, 1], label='Суммарная доля')
        
        plt.tight_layout()
        return fig
    
    def create_animation(self, save_path="virus_evolution.gif"):
        """
        Создает анимацию эволюции вирусной волны
        
        Parameters:
        -----------
        save_path : str
            Путь для сохранения анимации (gif)
        """
        if len(self.I_files) == 0:
            print("Нет данных для анимации!")
            return
            
        # Создаем фигуру
        self.fig, self.axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # Загружаем первый кадр для настройки цветовых шкал
        I0 = self.load_matrix(self.I_files[0])
        S0 = self.load_matrix(self.S_files[0])
        R0 = self.load_matrix(self.R_files[0])
        
        # Настраиваем subplots
        titles = ['Инфицированные (I)', 'Восприимчивые (S)', 'Выздоровевшие (R)']
        cmaps = ['hot', 'Blues', 'Greens']
        data = [I0, S0, R0]
        
        ims = []
        for i in range(3):
            im = self.axes[i].imshow(data[i], cmap=cmaps[i], aspect='auto', origin='lower', animated=True)
            self.axes[i].set_title(titles[i])
            self.axes[i].set_xlabel('Антигенная координата (x)')
            self.axes[i].set_ylabel('Нейтральная координата (y)')
            plt.colorbar(im, ax=self.axes[i])
            ims.append(im)
            
        self.fig.suptitle('Двумерная модель антигенной эволюции вируса\nШаг: 0', fontsize=14)
        plt.tight_layout()
        
        # Функция обновления кадра
        def update(frame):
            # Загружаем данные для текущего кадра
            I = self.load_matrix(self.I_files[frame])
            S = self.load_matrix(self.S_files[frame])
            R = self.load_matrix(self.R_files[frame])
            
            # Обновляем изображения
            ims[0].set_data(I)
            ims[1].set_data(S)
            ims[2].set_data(R)
            
            # Обновляем заголовок
            step_num = self.extract_step_number(self.I_files[frame])
            self.fig.suptitle(f'Двумерная модель антигенной эволюции вируса\nШаг: {step_num}', fontsize=14)
            
            # Обновляем цветовые шкалы
            for im in ims:
                im.autoscale()
                
            return ims
        
        # Создаем анимацию
        self.anim = FuncAnimation(self.fig, update, frames=len(self.I_files), interval=500, blit=True)
        
        # Сохраняем анимацию
        print(f"Создание анимации...")
        self.anim.save(save_path, writer=PillowWriter(fps=2))
        print(f"Анимация сохранена как {save_path}")
        
        return self.anim
    
    def plot_time_series(self):
        """Строит графики временных рядов нормировки и доли инфицированных"""
        norm_data, finf_data = self.load_time_series()
        
        if norm_data is None or finf_data is None:
            print("Файлы временных рядов не найдены!")
            return
            
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # График нормировки
        axes[0].plot(norm_data['Time'], norm_data['Norm'], 'b-', linewidth=2)
        axes[0].set_xlabel('Время')
        axes[0].set_ylabel('Нормировка (I+S)')
        axes[0].set_title('Нормировка популяции во времени')
        axes[0].grid(True, alpha=0.3)
        axes[0].axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Идеальное значение (1.0)')
        axes[0].legend()
        
        # График доли инфицированных
        axes[1].plot(finf_data['Time'], finf_data['Finf'], 'r-', linewidth=2)
        axes[1].set_xlabel('Время')
        axes[1].set_ylabel('Доля инфицированных')
        axes[1].set_title('Доля инфицированных во времени')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('output/time_series.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def plot_3d_surface(self, step_idx=-1):
        """
        Создает 3D поверхность для визуализации распределения инфицированных
        
        Parameters:
        -----------
        step_idx : int
            Индекс шага (по умолчанию -1 - последний шаг)
        """
        if step_idx == -1:
            step_idx = len(self.I_files) - 1
            
        if step_idx >= len(self.I_files):
            print(f"Ошибка: шаг {step_idx} не существует!")
            return
            
        # Загружаем данные
        I = self.load_matrix(self.I_files[step_idx])
        step_num = self.extract_step_number(self.I_files[step_idx])
        
        # Создаем 3D график
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Создаем координатную сетку
        x = np.arange(I.shape[1])
        y = np.arange(I.shape[0])
        X, Y = np.meshgrid(x, y)
        
        # Создаем поверхность
        surf = ax.plot_surface(X, Y, I, cmap=cm.hot, 
                              linewidth=0, antialiased=True, 
                              alpha=0.8, rstride=1, cstride=1)
        
        ax.set_xlabel('Антигенная координата (x)')
        ax.set_ylabel('Нейтральная координата (y)')
        ax.set_zlabel('Доля инфицированных')
        ax.set_title(f'3D распределение инфицированных\nШаг: {step_num}')
        
        # Добавляем цветовую шкалу
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Доля инфицированных')
        
        # Настраиваем угол обзора
        ax.view_init(elev=30, azim=45)
        
        plt.tight_layout()
        plt.savefig(f'output/3d_distribution_step_{step_num}.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def plot_contour_with_wave_front(self, step_idx=-1):
        """
        Создает контурный график с выделением фронта волны
        
        Parameters:
        -----------
        step_idx : int
            Индекс шага (по умолчанию -1 - последний шаг)
        """
        if step_idx == -1:
            step_idx = len(self.I_files) - 1
            
        if step_idx >= len(self.I_files):
            print(f"Ошибка: шаг {step_idx} не существует!")
            return
            
        # Загружаем данные
        I = self.load_matrix(self.I_files[step_idx])
        step_num = self.extract_step_number(self.I_files[step_idx])
        
        # Создаем фигуру
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # 1. Контурный график
        contour_levels = 20
        contour = axes[0].contourf(I, levels=contour_levels, cmap='hot')
        axes[0].set_title(f'Контурное распределение инфицированных\nШаг: {step_num}')
        axes[0].set_xlabel('Антигенная координата (x)')
        axes[0].set_ylabel('Нейтральная координата (y)')
        plt.colorbar(contour, ax=axes[0], label='Доля инфицированных')
        
        # 2. График с выделенным фронтом волны
        # Определяем порог для фронта волны (например, 50% от максимума)
        threshold = 0.5 * np.max(I)
        
        # Создаем маску для фронта волны
        wave_front = np.zeros_like(I)
        wave_front[I > threshold] = 1
        
        im = axes[1].imshow(I, cmap='gray', aspect='auto', origin='lower', alpha=0.7)
        axes[1].contour(wave_front, levels=[0.5], colors='red', linewidths=2)
        axes[1].set_title(f'Фронт волны инфицирования\nШаг: {step_num}')
        axes[1].set_xlabel('Антигенная координата (x)')
        axes[1].set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im, ax=axes[1], label='Доля инфицированных')
        
        plt.tight_layout()
        plt.savefig(f'output/wave_front_step_{step_num}.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def plot_all_steps_grid(self, max_steps=6):
        """
        Создает сетку графиков для нескольких шагов
        
        Parameters:
        -----------
        max_steps : int
            Максимальное количество шагов для отображения
        """
        n_steps = min(len(self.I_files), max_steps)
        
        # Выбираем равномерно распределенные шаги
        step_indices = np.linspace(0, len(self.I_files)-1, n_steps, dtype=int)
        
        fig, axes = plt.subplots(n_steps, 3, figsize=(15, 4*n_steps))
        
        # Если только один шаг, axes становится 1D массивом
        if n_steps == 1:
            axes = axes.reshape(1, -1)
        
        for i, step_idx in enumerate(step_indices):
            # Загружаем данные
            I = self.load_matrix(self.I_files[step_idx])
            S = self.load_matrix(self.S_files[step_idx])
            R = self.load_matrix(self.R_files[step_idx])
            
            step_num = self.extract_step_number(self.I_files[step_idx])
            
            # Отображаем матрицы
            im1 = axes[i, 0].imshow(I, cmap='hot', aspect='auto', origin='lower')
            axes[i, 0].set_title(f'Шаг {step_num}: Инфицированные')
            axes[i, 0].set_ylabel('y')
            if i == n_steps-1:
                axes[i, 0].set_xlabel('x')
            
            im2 = axes[i, 1].imshow(S, cmap='Blues', aspect='auto', origin='lower')
            axes[i, 1].set_title(f'Шаг {step_num}: Восприимчивые')
            if i == n_steps-1:
                axes[i, 1].set_xlabel('x')
            
            im3 = axes[i, 2].imshow(R, cmap='Greens', aspect='auto', origin='lower')
            axes[i, 2].set_title(f'Шаг {step_num}: Выздоровевшие')
            if i == n_steps-1:
                axes[i, 2].set_xlabel('x')
        
        plt.tight_layout()
        plt.savefig('output/all_steps_grid.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def analyze_wave_dynamics(self):
        """Анализирует динамику волны: скорость, ширину, форму"""
        if len(self.I_files) < 2:
            print("Недостаточно данных для анализа динамики!")
            return
            
        # Собираем данные о положении центра масс волны
        centers_x = []
        centers_y = []
        widths = []
        steps = []
        
        for i, I_file in enumerate(self.I_files):
            I = self.load_matrix(I_file)
            step_num = self.extract_step_number(I_file)
            
            if np.sum(I) > 0:
                # Вычисляем координаты центра масс
                y_coords, x_coords = np.indices(I.shape)
                center_x = np.sum(x_coords * I) / np.sum(I)
                center_y = np.sum(y_coords * I) / np.sum(I)
                
                # Вычисляем эффективную ширину
                variance_x = np.sum((x_coords - center_x)**2 * I) / np.sum(I)
                width = np.sqrt(variance_x)
                
                centers_x.append(center_x)
                centers_y.append(center_y)
                widths.append(width)
                steps.append(step_num)
        
        # Строим графики динамики
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # 1. Движение центра по x (антигенная координата)
        axes[0].plot(steps, centers_x, 'bo-', linewidth=2)
        axes[0].set_xlabel('Шаг времени')
        axes[0].set_ylabel('Центр волны по X')
        axes[0].set_title('Движение волны в антигенном направлении')
        axes[0].grid(True, alpha=0.3)
        
        # Линейная регрессия для оценки скорости
        if len(steps) > 1:
            coeffs = np.polyfit(steps, centers_x, 1)
            poly = np.poly1d(coeffs)
            axes[0].plot(steps, poly(steps), 'r--', 
                        label=f'Скорость: {coeffs[0]:.4f} клеток/шаг')
            axes[0].legend()
        
        # 2. Движение центра по y (нейтральная координата)
        axes[1].plot(steps, centers_y, 'go-', linewidth=2)
        axes[1].set_xlabel('Шаг времени')
        axes[1].set_ylabel('Центр волны по Y')
        axes[1].set_title('Движение волны в нейтральном направлении')
        axes[1].grid(True, alpha=0.3)
        
        # 3. Ширина волны
        axes[2].plot(steps, widths, 'ro-', linewidth=2)
        axes[2].set_xlabel('Шаг времени')
        axes[2].set_ylabel('Ширина волны')
        axes[2].set_title('Изменение ширины волны')
        axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('output/wave_dynamics.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        # Сохраняем результаты анализа
        analysis_data = pd.DataFrame({
            'Step': steps,
            'Center_X': centers_x,
            'Center_Y': centers_y,
            'Width': widths
        })
        analysis_data.to_csv('output/wave_analysis.csv', index=False)
        
        print("Анализ динамики волны завершен!")
        print(f"Средняя скорость по X: {coeffs[0]:.4f} клеток/шаг" if len(steps) > 1 else "")
        
        return fig, analysis_data

# Основная функция для запуска визуализации
def main():
    print("=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ДВУМЕРНОЙ МОДЕЛИ АНТИГЕННОЙ ЭВОЛЮЦИИ ВИРУСА")
    print("=" * 60)
    
    # Создаем визуализатор
    visualizer = VirusModelVisualizer("output")
    
    if len(visualizer.I_files) == 0:
        print("Нет данных для визуализации!")
        return
    
    # 1. Создаем сетку графиков для всех шагов
    print("\n1. Создание сетки графиков для всех шагов...")
    visualizer.plot_all_steps_grid(max_steps=6)
    
    # 2. Создаем анимацию
    print("\n2. Создание анимации эволюции...")
    try:
        visualizer.create_animation("virus_evolution.gif")
    except Exception as e:
        print(f"Ошибка при создании анимации: {e}")
    
    # 3. Строим графики временных рядов
    print("\n3. Построение графиков временных рядов...")
    visualizer.plot_time_series()
    
    # 4. Анализируем динамику волны
    print("\n4. Анализ динамики волны...")
    visualizer.analyze_wave_dynamics()
    
    # 5. Создаем 3D визуализацию для последнего шага
    print("\n5. Создание 3D визуализации...")
    visualizer.plot_3d_surface(step_idx=-1)
    
    # 6. Создаем контурные графики с фронтом волны
    print("\n6. Создание контурных графиков с фронтом волны...")
    visualizer.plot_contour_with_wave_front(step_idx=-1)
    
    # 7. Создаем отдельные графики для каждого шага
    print("\n7. Создание отдельных графиков для каждого шага...")
    for i in range(min(3, len(visualizer.I_files))):  # Первые 3 шага
        fig = visualizer.create_composite_plot(i)
        step_num = visualizer.extract_step_number(visualizer.I_files[i])
        plt.savefig(f'output/composite_plot_step_{step_num}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
    
    print("\n" + "=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print("Все графики сохранены в папке 'output'")
    print("=" * 60)
    
    # Выводим сводную информацию
    print("\nСводная информация:")
    print(f"  • Всего шагов симуляции: {len(visualizer.I_files)}")
    print(f"  • Размер решётки: {visualizer.load_matrix(visualizer.I_files[0]).shape}")
    print(f"  • Созданные файлы:")
    print(f"    - Анимация: virus_evolution.gif")
    print(f"    - Временные ряды: output/time_series.png")
    print(f"    - Динамика волны: output/wave_dynamics.png")
    print(f"    - 3D визуализация: output/3d_distribution_step_*.png")
    print(f"    - Контурные графики: output/wave_front_step_*.png")
    print(f"    - Сетка графиков: output/all_steps_grid.png")

if __name__ == "__main__":
    # Создаем папку для графиков, если её нет
    if not os.path.exists("output"):
        os.makedirs("output")
    
    main()