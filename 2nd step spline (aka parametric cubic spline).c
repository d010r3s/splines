#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*На оценку 3 достаточно реализовать построение одного сплайна по
заданным точкам. Иными словами, нужно, имея массив из некоторо-
го количества точек, написать программу, которая позволит определять
значение функции в любой промежуточной точке. При этом гарантиру-
ется, что кривая, проходящая через данный массив точек, является од-
нозначной функцией, то есть каждому значению x соответствует только
одно значение y.
1 балл добавляется, если сплайн построен параметрическим образом,
то есть допускает задание неоднозначных функций (произвольных кри-
вых, когда одному значению x соответствует несколько значений y).
Плюс 1 балл дается за реализацию алгоритма нахождения точки пе-
ресечения двух сплайнов. Еще 1 балл добавляется, если реализованный
алгоритм признается достойным студента программы ПМИ (в достаточ-
ной мере использует материалы математических курсов).
Плюс 1 балл дается за реализацию алгоритма нахождения наимень-
шего расстояния между двумя сплайнами. */


double spline(double x, double * node, double * val, int n)
{
    int i;
    double h[n], b[n], u[n], v[n], z[n], c[n], d[n], p;
    // массив расстояний между узлами
    for (i = 0; i < n - 1; i++) {
        h[i] = node[i+1] - node[i];
    }

    // массив расстояний между значениями
    for (i = 0; i < n - 1; i++) {
        b[i] = (val[i+1] - val[i]) / h[i];
    }

    // метод прогонки (можно почитать тут: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
    u[1] = 2.0 * (h[0] + h[1]);
    v[1] = 6.0 * (b[1] - b[0]);
    for (i = 2; i < n - 1; i++) {
        u[i] = 2.0 * (h[i] + h[i-1]) - (h[i-1] * h[i-1]) / u[i-1];
        v[i] = 6.0 * (b[i] - b[i-1]) - (h[i-1] * v[i-1]) / u[i-1];
    }
    z[n-1] = 0.0;
    for (i = n - 2; i > 0; i--) {
        z[i] = (v[i] - h[i] * z[i+1]) / u[i];
    }
    z[0] = 0.0;

    // коэффициенты сплайна
    for (i = 0; i < n - 1; i++) {
        c[i] = (val[i+1] - val[i]) / h[i] - h[i] * (z[i+1] + 2.0 * z[i]) / 6.0;
        d[i] = z[i] / 2.0;
    }

    // значение сплайна в точке
    //https://www.physicsforums.com/attachments/parametric-spline-tutorialv2-pdf.12898/ почитать тут
    for (i = 0; i < n - 1; i++) {
        if (x >= node[i] && x <= node[i+1]) {
            p = c[i] + (x - node[i]) * (d[i] + (x - node[i]) * (z[i] - d[i]) / h[i]) / h[i];
            return p;
        }
    }
}

int main()
{
    int n;
    printf("Enter the number of knots:");
    scanf("%d", &n);
    double x;
    double * node = malloc(sizeof(double)*n);
    printf("Enter the knots:");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &node[i]);
    }
    double * val = malloc(sizeof(double) * n);
    printf("Enter the values:");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &val[i]);
    }
    printf("Enter the value of x:");
    scanf("%lf", &x);

    double result = spline(x, node, val, n);

    printf("The value of the parametric cubic spline at x = %lf is %lf\n", x, result);
}
