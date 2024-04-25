import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('output.txt', usecols=(0, 1, 2), max_rows=100000)
x = data[:, 0]
y = data[:, 1]
t = data[:, 2]

plt.figure()

plt.subplot(2, 1, 1)
plt.plot(t, x, 'r', label='x(t)')
plt.plot(t, y, 'b', label='y(t)')
plt.legend()
plt.xlabel('t')
plt.ylabel('Значение')

plt.subplot(2, 1, 2)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')

plt.show()