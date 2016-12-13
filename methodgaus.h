#ifndef METHODGAUS_H
#define METHODGAUS_H

#include <QVector>

class methodGaus
{

private:
    QVector<QVector<double>> matrix;
    QVector<double> answers;
    const double eps = 10e-12;

public:
    methodGaus(QVector<QVector<double>> matrix) : matrix(matrix)
    {
        answers.resize(matrix.size());
    }

public:
   QVector<double> calculate() throw (QString);

private:
   int  searchMaxInCurrentColumn(int currentRow);

};

#endif // METHODGAUS_H
