#include "mainwindow.h"
#include "ui_mainwindow.h"

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
        qDebug() << counter;
        if (counter % Resources::MAX_ITERS == 0) {
            break;
        }
    }
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
