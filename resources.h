#include <qmath.h>
#include <interpolation.h>

typedef QVector<int> Ints;
typedef QVector<double> Doubles;
typedef QVector<QVector<double>> Doubles2D;

V_Nodes absorb_nodes = {
    QPointF(293, 2e-3),
    QPointF(1278, 5e-3),
    QPointF(1528, 7.8e-3),
    QPointF(1677, 10e-3),
};

struct Resources {
    // Теплоёмкость
    static double calc_C(double T) {
        return 2.049 + 0.563e-3 * T - 0.528e5 / T / T;
    }

    // Коэффициент теплопроводности
    static double calc_lambda(double T) {
        return 1.34e-1 * (1 + 4.35e-4 * T);
    }

    // Коэффициент преломления
    constexpr static double REFRACTION_F = 1.46;

    // Коэффициент поглощения - интерполируется по таблице
    static double calc_absorb(double T) throw (QString) {
        return Interpolation().forward(T, 2, absorb_nodes);
    }

    // Функция f из условия
    static double calc_f(double F0, double Z, double T) throw (QString) {
        return F0 * exp(-calc_absorb(T) * Z);
    }

    static double init_lambdas(const double U0, Doubles2D & lambdas) throw (QString) {
        double start_lambda = Resources::calc_lambda(U0);
        for (Doubles & lambda_row : lambdas) {
            for (int J = 0; J < lambda_row.size(); J++) {
                lambda_row[J] = start_lambda;
            }
        }
    }

};
