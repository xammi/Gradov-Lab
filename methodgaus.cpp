#include "methodgaus.h"
#include <qmath.h>

QVector<double> methodGaus::calculate()
    throw (QString)
{
    int currentRow = 0;

    while(currentRow < matrix.size())
    {
        int maxIndexRow = searchMaxInCurrentColumn(currentRow);

        //Вырожденный случай, что данная переменная может принимать произвольное значение
        //-----------------------------------------------
        if (fabs (matrix[maxIndexRow][currentRow]) < eps)
             throw QString("Решение получить невозможно из-за нулевого диагонального элемента");
        //-----------------------------------------------
        //Меняем местами строки в матрице
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
            answers[i] = answers[i] - matrix[i][currentRow] * answers[currentRow];
    }

    return answers;
}

int methodGaus::searchMaxInCurrentColumn(int currentRow)
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
