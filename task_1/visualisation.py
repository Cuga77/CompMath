import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Размер окна
WINDOW_WIDTH = 125
WINDOW_HEIGHT = 125

# Чтение данных из файла
data = np.loadtxt("output.txt")

# Разделение данных на x, y, vx, vy
x = data[:, 0].reshape((-1, WINDOW_WIDTH, WINDOW_HEIGHT))
y = data[:, 1].reshape((-1, WINDOW_WIDTH, WINDOW_HEIGHT))
vx = data[:, 2].reshape((-1, WINDOW_WIDTH, WINDOW_HEIGHT))
vy = data[:, 3].reshape((-1, WINDOW_WIDTH, WINDOW_HEIGHT))


fig = plt.figure()

# Функция для обновления изображения на каждом шаге
def update(i):
    plt.clf()
    plt.imshow(x[i], cmap='winter_r', interpolation='none')
    # plt.clf()
    # plt.imshow(y[i], cmap='winter', interpolation='none')
    # plt.clf()
    # plt.imshow(vy[i], cmap='winter', interpolation='none')
    # plt.clf()
    # plt.imshow(vx[i], cmap='winter', interpolation='none')

ani = animation.FuncAnimation(fig, update, frames=range(len(x)), interval=1)

plt.show()