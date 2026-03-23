import numpy as np
import matplotlib.pyplot as plt
import os

stept = 0.5
Tmax = 2000.0
L = 120
data_dir = "out"

# Загружаем нормировку и долю инфицированных
norm = np.loadtxt(os.path.join(data_dir, "norm.txt"))
finf = np.loadtxt(os.path.join(data_dir, "finf.txt"))
time = np.arange(stept, Tmax + 1e-9, stept)

# Загружаем моменты времени сохранённых срезов
slice_times = np.loadtxt(os.path.join(data_dir, "slice_times.txt"))
n_slices = len(slice_times)

# Подготовка фигуры
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Левый график: все срезы I и R
x_vals = np.arange(L)

for idx, t in enumerate(slice_times):
    I_slice = np.loadtxt(os.path.join(data_dir, f"I_slice_{idx}.txt"))
    R_slice = np.loadtxt(os.path.join(data_dir, f"R_slice_{idx}.txt"))
    
    # Прозрачность меняется от 0.3 до 1.0, чтобы различать слои
    alpha = 0.3 + 0.7 * idx / (n_slices - 1) if n_slices > 1 else 1.0
    
    # Инфекция (I) – красный, восстановленные (R) – синий
    ax1.plot(x_vals, 10 * I_slice, '-', color='red', alpha=alpha, label='10*I' if idx == 0 else "")
    ax1.plot(x_vals, R_slice, '-', color='blue', alpha=alpha, label='R' if idx == 0 else "")

ax1.set_xlabel('Mutation number')
ax1.set_ylabel('Density')
ax1.set_title('Distributions at different times')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend()

# Правый график: доля инфицированных во времени
ax2.plot(time, finf, 'k-')
ax2.set_xlabel('Time')
ax2.set_ylabel('Fraction infected (finf)')
ax2.set_title('Fraction infected over time')
ax2.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.show()