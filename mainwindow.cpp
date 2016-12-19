#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <time.h>

typedef struct timespec Time;

Time get_current_time() {
    Time gettime_now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, & gettime_now);
    return gettime_now;
}

void get_end_time(Time start_time) {
    Time end_time = get_current_time();
    double sec = double(end_time.tv_sec - start_time.tv_sec);
    double nsec = double(end_time.tv_nsec - start_time.tv_nsec) / 10e9;
    sec -= nsec;

    int hours = qFloor(sec) / 3600;
    int minutes = (qFloor(sec) % 3600) / 60;
    int seconds = (qFloor(sec) % 3600) % 60;
    int millis = qFloor(sec * 1000) % 1000;

    QString timing = "Затраченное время:";
    if (hours > 0) timing += QString::number(hours) + " часов ";
    if (minutes > 0) timing += QString::number(minutes) + " минут ";
    if (seconds > 0) timing += QString::number(seconds) + "." + QString::number(millis) + " секунд";

    qDebug() << timing;
}
//-------------------------------------------------------------------------------------------------

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    gauss(new GaussResolver())
{
    ui->setupUi(this);
    connect(ui->btn_recalc, SIGNAL(pressed()), SLOT(recalculate()));
}

MainWindow::~MainWindow() {
    delete ui;
    delete gauss;
}

void MainWindow::recalculate() {
    ui->T->setText("");
    ui->error->setText("");
    try {
        Doubles2D Ts;
        this->straight_method(Ts);
        this->view_result(Ts);
    }
    catch (QString error) {
        ui->error->setText(error);
    }
}

void MainWindow::straight_method(Doubles2D & Ts) throw (QString) {
    double A = ui->A->value(), B = ui->B->value();
    int Nx = ui->Nx->value(), Nz = ui->Nz->value();
    double Ft = ui->Ft->value(), F0 = ui->F0->value();
    double U0 = ui->U0->value(), alpha = ui->alpha->value();

    if (Nx == 0) {
        throw "Количество узлов сетки по X равно 0";
    }
    if (Nz == 0) {
        throw "Количество узлов сетки по Z равно 0";
    }
    gauss->initialize(Nx, Nz);

    double Hx = A / Nx, Hz = B / Nz;
    double Hx2 = Hx * Hx, Hz2 = Hz * Hz;
    int dim = Nx * Nz;
    int Nx_1 = Nx - 1, Nz_1 = Nz - 1;

    Doubles2D lambdas;
    rs.init_lambdas(U0, Nx, Nz, lambdas);
    Doubles2D prev_Ts;
    rs.init_Ts(Nx, Nz, Ts);
    Doubles2D matrix;
    rs.init_matrix(Nx, Nz, matrix);

    Time start_time = get_current_time();
    int counter = 0;
    while (rs.if_stop_iterations(Nx, Nz, Ts, prev_Ts)) {
        // вычисление 1-го краевого условия (T[0][J] - T[1][J] = K1[J])
        for (int J = 0; J < Nx; J++) {
            matrix[J][dim] = Hz * Ft / lambdas[0][J];
        }

        // вычисление 2-го краевого условия (K1[I] * T[I][0] - T[I][1] = K2[I])
        for (int I = 1; I < Nz; I++) {
            matrix[I*Nx][dim] = -alpha * Hx / lambdas[I][0] * U0;
        }

        // вычисление 3-го краевого условия (T[I][Nx - 1] - K1[I] * T[I][Nx] = K2[I])
        for (int I = 1; I < Nz_1; I++) {
            matrix[(I + 1)*Nx - 1][dim] = -alpha * Hx / lambdas[I][Nx_1] * U0;
        }

        // вычисление 4-го краевого условия (T[Nz - 1][J] - K1[J] * T[Nz][J] = K2[J])
        for (int J = 1; J < Nx; J++) {
            matrix[Nz_1*Nx + J][dim] = -alpha * Hz / lambdas[Nz_1][J] * U0;
        }

        // заполнение правой части
        for (int I = 1; I < Nz_1; I++) {
            for (int J = 1; J < Nx_1; J++) {
                matrix[I*Nx + J][dim] = -rs.calc_f(F0, Hz * I, Ts[I][J]);
            }
        }

        for (int J = 0; J < Nx; J++) {
            matrix[J][J] = 1;
            matrix[J][Nx + J] = -1;
        }
        for (int I = 1; I < Nz_1; I++) {
            matrix[I*Nx][I*Nx] = 1 + alpha * Hx / lambdas[I][0];
            matrix[I*Nx][I*Nx + 1] = -1;
        }
        for (int I = 1; I < Nz_1; I++) {
            matrix[I*Nx + Nx_1][I*Nx + Nx_1 - 1] = 1;
            matrix[I*Nx + Nx_1][I*Nx + Nx_1] = -(1 + alpha * Hx / lambdas[I][Nx_1]);
        }
        for (int J = 0; J < Nx; J++) {
            matrix[Nz_1*Nx + J][(Nz_1 - 1)*Nx + J] = 1;
            matrix[Nz_1*Nx + J][Nz_1*Nx + J] = -(1 + alpha * Hz / lambdas[Nz_1][J]);
        }

        // заполнение диагоналей в центре
        for (int I = 1; I < Nz_1; I++) {
            for (int J = 1; J < Nx_1; J++) {
                double bottom_L = (lambdas[I][J] + lambdas[I][J - 1]) / 2;
                double top_L = (lambdas[I][J] + lambdas[I][J + 1]) / 2;
                double left_L = (lambdas[I][J] + lambdas[I - 1][J]) / 2;
                double right_L = (lambdas[I][J] + lambdas[I + 1][J]) / 2;

                matrix[I*Nx + J][(I - 1)*Nx + J] = left_L / Hz2;
                matrix[I*Nx + J][I*Nx + (J - 1)] = bottom_L / Hx2;

                matrix[I*Nx + J][I*Nx + J] = -((top_L + bottom_L) / Hx2 + (left_L + right_L) / Hz2);

                matrix[I*Nx + J][I*Nx + (J + 1)] = top_L / Hx2;
                matrix[I*Nx + J][(I + 1)*Nx + J] = right_L / Hz2;
            }
        }

        prev_Ts = Ts;
        gauss->set_matrix(matrix);
        Doubles & solution = gauss->calculate();
        for (int I = 0; I < Nz; I++) {
            for (int J = 0; J < Nx; J++) {
                Ts[I][J] = solution[I * Nz + J];
            }
        }
        rs.recalc_lambdas(Nx, Nz, Ts, lambdas);

        counter++;
        qDebug() << "Итерация №:" << counter;
        if (counter % Resources::MAX_ITERS == 0) {
            break;
        }
    }
    get_end_time(start_time);
}

void MainWindow::relaxation_method(Doubles2D & Ts) throw (QString) {
    double A = ui->A->value(), B = ui->B->value();
    int Nx = ui->Nx->value(), Nz = ui->Nz->value();
    double Ft = ui->Ft->value(), F0 = ui->F0->value();
    double U0 = ui->U0->value(), alpha = ui->alpha->value();
    double tau = 0.01;

    if (Nx == 0) {
        throw "Количество узлов сетки по X равно 0";
    }
    if (Nz == 0) {
        throw "Количество узлов сетки по Z равно 0";
    }

    double Hx = A / Nx, Hz = B / Nz;
    double Hx2 = Hx * Hx, Hz2 = Hz * Hz;
    int Nx_1 = Nx - 1, Nz_1 = Nz - 1, Nx_2 = Nx - 2, Nz_2 = Nz - 2;

    Time start_time = get_current_time();

    Doubles2D lambdas;
    rs.init_lambdas(U0, Nx, Nz, lambdas);
    Doubles2D Ts_strich, Ts_roof, prev_Ts;
    rs.init_Ts(Nx, Nz, Ts); rs.init_Ts(Nx, Nz, Ts_strich); rs.init_Ts(Nx, Nz, Ts_roof);

    Doubles alphas_I, betas_I;
    alphas_I.fill(0.0, Nx); betas_I.fill(0.0, Nx);
    auto A1 = [&lambdas, Hx2] (I, J) {
        return (lambdas[I][J] + lambdas[I][J + 1]) / 2 / Hx2;
    };
    auto C1 = [&lambdas, Hx2, tau] (I, J) {
        return (-1) * (2 * lambdas[I][J] + lambdas[I][J + 1] + lambdas[I][J - 1]) / 2 / Hx2 + 1 / tau;
    };
    auto B1 = [&lambdas, Hx2] (I, J) {
        return (lambdas[I][J] + lambdas[I][J - 1]) / 2 / Hx2;
    };
    auto F1 = [this, &Ts, F0, tau, Hz] (I, J) {
        return (-1) * this->rs.calc_f(F0, Hz * I, Ts[I][J]) / 2 + Ts[I][J] / tau;
    };

    Doubles alphas_J, betas_J;
    alphas_J.fill(0.0, Nz); betas_J.fill(0.0, Nz);

    int counter = 0;
    while (rs.if_stop_iterations(Nx, Nz, Ts, prev_Ts)) {
        prev_Ts = Ts;
        if (counter > 0) {
            Ts = Ts_roof;
            rs.recalc_lambdas(Nx, Nz, Ts, lambdas);
        }

        // прогонка для фиксированных I
        for (int I = 1; I < Nz_1; I++) {
            double C_I1 = C1(I, 1);
            alphas_I[2] = -B1(I, 1) / C_I1;
            betas_I[2] = F1(I, 1) / C_I1;

            // вычисление прогоночных коэффициентов
            for (int J = 1; J < Nx_1; J++) {
                double A_J = A1(I, J);
                double znam = (A_J * alphas_I[J] + C1(I, J));
                alphas_I[J + 1] = -B1(I, J) / znam;
                betas_I[J + 1] = (F1(I, J) - A_J * betas_I[J]) / znam;
            }
            double A_Nx = A1(I, Nx_2);
            Ts_strich[I][Nx_2] = (F1(I, Nx_2) - A_Nx * betas_I[Nx_2]) / (C1(I, Nx_2) + A_Nx * alphas_I[Nx_2]);

            // по краевому условию (4)
            double K4 = alpha * Hz / lambdas[I][Nx_1];
            Ts_strich[I][Nx_1] = Ts_strich[I][Nx_2] + K4 * U0 / (1 + K4);

            // восстановление Y с надчёркиванием
            for (int J = Nx_2 - 1; J > 0; J--) {
                Ts_strich[I][J] = alphas_I[J] * Ts_strich[I][J + 1] + betas_I[J];
            }

            // по краевому условию (3)
            double K3 = alpha * Hx / lambdas[I][0];
            Ts_strich[I][0] = (Ts_strich[I][1] - K3 * U0) / (1 + K3);
        }
        for (int J = 0; J < Nx; J++) {
            // по краевому условию (2)
            Ts_strich[0][J] = Ts_strich[1][J] + Hz / lambdas[0][J] * Ft;

            // по краевому условию (5)
            double K5 = alpha * Hz / lambdas[Nz_1][J];
            Ts_strich[Nz_1][J] = Ts_strich[Nz_2][J] * (1 + K5) - K5 * U0;
        }

        // resolve to Ts_roof basing on Ts_strich

        counter++;
        qDebug() << "Итерация №:" << counter;
        if (counter % Resources::MAX_ITERS == 0) {
            break;
        }
    }
    get_end_time(start_time);
}

void MainWindow::view_result(Doubles2D &Ts) {
    double A = ui->A->value(), B = ui->B->value();
    int Nx = ui->Nx->value(), Nz = ui->Nz->value();
    double Hx = A / Nx, Hz = B / Nz;

    ui->T->append("\"x\", \"z\", \"T\"");

    int G = 0, I = 0, J = 0;
    for (double Z = 0.0; Z < A; Z += Hz) {
        J = 0;
        for (double X = 0.0; X < B; X += Hx) {
            if (I < Nz && J < Nx) {
                ui->T->append(QString::number(Z) + ", " +
                              QString::number(X) + ", " +
                              QString::number(Ts[I][J]));
            }
            G++; J++;
        }
        I++;
    }
}
