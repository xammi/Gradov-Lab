#ifndef METHODGAUS_H
#define METHODGAUS_H

#include <QVector>

class GaussResolver
{
private:
    const double eps = 1e-80;

    QVector<QVector<double>> matrix;
    QVector<double> answers;
    int dim;

public:
    void initialize(int Nx, int Nz);
    void set_matrix(QVector<QVector<double>> & data);

    QVector<double> & calculate_v1() throw (QString);
    QVector<double> & calculate();

private:
    int  search_max_in_column(int currentRow);
    void forward_gauss();
    void backward_gauss();

};

#endif // METHODGAUS_H
