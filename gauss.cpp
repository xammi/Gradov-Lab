#include "gauss.h"
#include <qmath.h>

QVector<double> & GaussResolver::calculate_v1()
    throw (QString)
{
    int currentRow = 0;

    while(currentRow < dim)
    {
        int maxIndexRow = search_max_in_column(currentRow);

        // Вырожденный случай, что данная переменная может принимать произвольное значение
        //-----------------------------------------------
        if (fabs (matrix[maxIndexRow][currentRow]) < eps)
             throw QString("Решение получить невозможно из-за нулевого диагонального элемента");
        //-----------------------------------------------
        // Меняем местами строки в матрице
        matrix[maxIndexRow].swap(matrix[currentRow]);

         //----Нормализация
        for (int i = currentRow; i < matrix.size(); i++)
        {
            double temp = matrix[i][currentRow];
            if (fabs(temp) < eps)
                continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < matrix[0].size(); j++)
            {
                matrix[i][j] = matrix[i][j] / temp;
            }
            if (i == currentRow)
                continue; // уравнение не вычитать само из себя
            for (int j = 0; j < matrix[0].size(); j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[currentRow][j];
            }
         }

        currentRow++;
    }

//------Обратный ход
    for(currentRow = matrix.size() - 1; currentRow >= 0; currentRow--)
    {
        answers[currentRow] = matrix[currentRow][matrix[0].size()-1];

        for(int i = 0; i < currentRow; i++)
            matrix[i][matrix[0].size()-1] = matrix[i][matrix[0].size()-1] - matrix[i][currentRow] * answers[currentRow];
    }

    return answers;
}

int GaussResolver::search_max_in_column(int currentRow)
{
    int maxIndexRow = currentRow;
    double maxElement = fabs(matrix[currentRow][currentRow]);

    for (int i = currentRow + 1; i < matrix.size(); i++)
        if (fabs(matrix[i][currentRow]) > maxElement)
        {
            maxElement = fabs(matrix[i][currentRow]);
            maxIndexRow = i;
        }

    return maxIndexRow;
}
//-------------------------------------------------------------------------------------------------
void GaussResolver::initialize(int Nx, int Nz) {
    dim = Nx * Nz;
    double start_val = 0.0;
    QVector<double> matrix_row;
    matrix_row.fill(start_val, dim + 1);
    matrix.fill(matrix_row, dim);
    answers.fill(start_val, dim);
}

void GaussResolver::set_matrix(QVector<QVector<double>> & data) {
    for (int I = 0; I < dim; I++) {
        for (int J = 0; J < dim + 1; J++) {
            matrix[I][J] = data[I][J];
        }
    }
}

QVector<double> & GaussResolver::calculate() {
    forward_gauss();
    backward_gauss();
    return answers;
}

void GaussResolver::forward_gauss() {
    double v;
    for (int k = 0,i,j,im; k < dim - 1; k++) {
        im = k;
        for (i = k + 1; i < dim; i++) {
            if (matrix[i][k] != 0) {
                if (fabs(matrix[im][k]) < fabs(matrix[i][k])) {
                    im = i;
                }
            }
        }
        if (im != k) {
            for (j = 0; j < dim; j++) {
                if (matrix[im][j] != matrix[k][j]) {
                    v = matrix[im][j];
                    matrix[im][j] = matrix[k][j];
                    matrix[k][j] = v;
                }
            }
            v = matrix[im][dim];
            matrix[im][dim] = matrix[k][dim];
            matrix[k][dim] = v;
        }
        for (i = k + 1; i < dim; i++) {
            v = 1.0 * matrix[i][k] / matrix[k][k];
            matrix[i][k] = 0;
            if (matrix[k][dim] != 0) {
                matrix[i][dim] = matrix[i][dim] - v * matrix[k][dim];
            }
            if (v != 0) {
                for (j = k + 1; j < dim; j++) {
                    if (matrix[k][j] != 0) {
                        matrix[i][j] = matrix[i][j] - v * matrix[k][j];
                    }
                }
            }
        }
    }
}

void GaussResolver::backward_gauss() {
    double s = 0;
    answers[dim - 1] = 1.0 * matrix[dim - 1][dim] / matrix[dim - 1][dim - 1];
    for (int i = dim - 2, j; 0 <= i; i--) {
        s = 0;
        for (j = i + 1; j < dim; j++) {
            if (matrix[i][j] != 0) {
                s = s + matrix[i][j] * answers[j];
            }
        }
        answers[i] = 1.0 * (matrix[i][dim] - s) / matrix[i][i];
    }
}
