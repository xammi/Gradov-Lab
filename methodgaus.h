#ifndef METHODGAUS_H
#define METHODGAUS_H

#include <QVector>

class methodGaus
{

private:
    const double eps = 1e-16;

    QVector<QVector<double>> matrix;
    QVector<double> answers;
    int dim;

public:
    void initialize(int Nx, int Nz) {
        dim = Nx * Nz;
        double start_val = 0.0;
        QVector<double> matrix_row;
        matrix_row.fill(start_val, dim + 1);
        matrix.fill(matrix_row, dim);
        answers.fill(start_val, dim);
    }

    void set_matrix(QVector<QVector<double>> & data)
    {
        for (int I = 0; I < dim; I++) {
            for (int J = 0; J < dim + 1; J++) {
                matrix[I][J] = data[I][J];
            }
        }
    }

public:
   QVector<double> & calculate() throw (QString);

private:
   int  searchMaxInCurrentColumn(int currentRow);

};

#endif // METHODGAUS_H
