#include <GLFW/glfw3.h>

int main() {
    // Инициализация библиотеки
    if (!glfwInit()) {
        return -1;
    }

    // Создание окна
    GLFWwindow* window = glfwCreateWindow(640, 480, "OpenGL", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Установка контекста окна
    glfwMakeContextCurrent(window);

    // Главный цикл
    while (!glfwWindowShouldClose(window)) {
        // Очистка буфера перед рисованием
        glClear(GL_COLOR_BUFFER_BIT);

        // Начало рисования
        glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f); // Установка цвета линии
        glVertex2f(-1.0f, -1.0f); // Начальная точка линии
        glVertex2f(1.0f, 1.0f); // Конечная точка линии
        glEnd();

        // Обмен буферов
        glfwSwapBuffers(window);

        // Проверка событий
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}