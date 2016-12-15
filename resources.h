#ifndef RESOURCES_H
#define RESOURCES_H

#include <qmath.h>
#include "interpolation.h"

typedef QVector<int> Ints;
typedef QVector<double> Doubles;
typedef QVector<QVector<double>> Doubles2D;


struct Resources {
    V_Nodes absorb_nodes;
    Resources() {
        absorb_nodes = {
            QPointF(293, 2e-3),
            QPointF(1278, 5e-3),
            QPointF(1528, 7.8e-3),
            QPointF(1677, 10e-3),
        };
    }

    // Теплоёмкость
    double calc_C(double T) {
        return 2.049 + 0.563e-3 * T - 0.528e5 / T / T;
    }

    // Коэффициент теплопроводности
    double calc_lambda(double T) {
        return 0.134e1 * (1 + 4.35e-4 * T);
    }

    // Коэффициент преломления
    constexpr static double REFRACTION_F = 1.46;

    // Коэффициент поглощения - интерполируется по таблице
    double calc_absorb(double T) throw (QString) {
        return Interpolation().forward(T, 4, absorb_nodes);
    }

    // Функция f из условия
    double calc_f(double F0, double Z, double T) throw (QString) {
        return F0 * exp(-calc_absorb(T) * Z);
    }

    void init_lambdas(double U0, int Nx, int Nz, Doubles2D & lambdas) throw (QString) {
        double start_lambda = this->calc_lambda(U0);
        Doubles lambda_row;
        lambda_row.fill(start_lambda, Nx);
        lambdas.fill(lambda_row, Nz);
    }

    void init_Ts(int Nx, int Nz, Doubles2D & Ts) {
        double start_T = 293.0;
        Doubles Ts_row;
        Ts_row.fill(start_T, Nx);
        Ts.fill(Ts_row, Nz);
    }

    void init_matrix(int Nx, int Nz, Doubles2D & matrix) {
        double start_val = 0.0;
        int dim = Nx * Nz;

        Doubles matrix_row;
        matrix_row.fill(start_val, dim + 1);
        matrix.fill(matrix_row, dim);
    }

    void recalc_lambdas(int Nx, int Nz, const Doubles2D & Ts, Doubles2D & lambdas) throw (QString) {
        for (int I = 0; I < Nz; I++) {
            for (int J = 0; J < Nx; J++) {
                lambdas[I][J] = this->calc_lambda(Ts[I][J]);
            }
        }
    }

    static constexpr double EPS = 10e-5;
    static constexpr int MAX_ITERS = 50;

    bool if_stop_iterations(int Nx, int Nz, Doubles2D & Ts, Doubles2D & prev_Ts) {
        if (prev_Ts.empty()) {
            return true;
        }
        for (int I = 0; I < Nz; I++) {
            for (int J = 0; J < Nx; J++) {
                if (fabs(prev_Ts[I][J] - Ts[I][J]) > Resources::EPS)
                    return true;
            }
        }
        return false;
    }
};

#endif // RESOURCES_H
