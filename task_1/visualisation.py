import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

WINDOW_WIDTH = 100
WINDOW_HEIGHT = 100

def animate(i):
    with open("./output.txt", "r") as file:
        lines = file.readlines()
        concentration = np.array([float(x) for x in lines[i].split()]).reshape((WINDOW_WIDTH, WINDOW_HEIGHT))

    # Визуализация
    im.set_array(concentration)
    return [im]

fig, ax = plt.subplots()
im = ax.imshow(np.ones((WINDOW_WIDTH, WINDOW_HEIGHT)), animated=True)
ani = animation.FuncAnimation(fig, animate, frames=100, interval=100, blit=True)

plt.show()