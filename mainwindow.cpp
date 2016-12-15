#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->btn_recalc, SIGNAL(pressed()), SLOT(recalculate()));
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::recalculate() {
    ui->T->setText("");
    ui->error->setText("");
    ui->pb->setValue(0);
    try {
        Doubles2D Ts;
        this->recalculate_action(Ts);
        ui->pb->setValue(100);
        this->view_result(Ts);
    }
    catch (QString error) {
        ui->error->setText(error);
    }
}

void MainWindow::recalculate_action(Doubles2D & Ts) throw (QString) {
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
        ui->pb->setValue(5);

        // вычисление 1-го краевого условия (T[0][J] - T[1][J] = K1[J])
        for (int J = 0; J < Nx; J++) {
            double K1 = Hz * Ft / lambdas[0][J];
            matrix[J][dim] = K1;
        }

        // вычисление 2-го краевого условия (K1[I] * T[I][0] - T[I][1] = K2[I])
        for (int I = 1; I < Nz - 1; I += 2) {
            double K1 = -alpha * Hx / lambdas[I][0] * U0;
            matrix[Nx * I][dim] = K1;
        }

        // вычисление 3-го краевого условия (T[I][Nx - 1] - K1[I] * T[I][Nx] = K2[I])
        for (int I = 1; I < Nz - 1; I += 2) {
            double K1 = -alpha * Hx / lambdas[I][Nx - 1] * U0;
            matrix[Nx + Nx * I][dim] = K1;
        }

        // вычисление 4-го краевого условия (T[Nz - 1][J] - K1[J] * T[Nz][J] = K2[J])
        for (int J = 0; J < Nx; J++) {
            double K1 = -alpha * Hz / lambdas[Nz - 1][J] * U0;
            matrix[Nz_1*Nx + J][dim] = K1;
        }
        ui->pb->setValue(20);

        // заполнение правой части
        for (int I = 1; I < Nz - 1; I++) {
            for (int J = 1; J < Nx - 1; J++) {
                matrix[I * Nx + J][dim] = -rs.calc_f(F0, Hz * I, Ts[I][J]);
            }
        }
        ui->pb->setValue(40);

        for (int J = 0; J < Nx; J++) {
            matrix[J][J] = 1;
            matrix[J][Nx + J] = -1;
            matrix[Nx + J][J] = -1;
        }

        // заполнение диагоналей в центре
        double augenblick;
        for (int I = 1; I < Nz - 1; I++) {
            for (int J = 0; J < Nx; J++) {
                if (J == 0) { augenblick = lambdas[I][J]; } else { augenblick = lambdas[I][J - 1]; }
                double bottom_L = (lambdas[I][J] + augenblick) / 2;

                if (J == Nx - 1) { augenblick = lambdas[I][J]; } else { augenblick = lambdas[I][J + 1]; }
                double top_L = (lambdas[I][J] + augenblick) / 2;

                double left_L = (lambdas[I][J] + lambdas[I - 1][J]) / 2;
                double right_L = (lambdas[I][J] + lambdas[I + 1][J]) / 2;

                if (J == 0) {
                    matrix[I*Nx + J][(I - 1)*Nx + J] = 1 - alpha * Hx / lambdas[I][0];
                    matrix[I*Nx + J][I*Nx + J - 1] = -1;
                } else {
                    matrix[I*Nx + J][(I - 1)*Nx + J] = left_L / Hz2;
                    matrix[I*Nx + J][I*Nx + J - 1] = bottom_L / Hx2;
                }

                matrix[I*Nx + J][I*Nx + J] = (top_L + bottom_L) / Hx2 + (left_L + right_L) / Hz2;

                if (J == Nx - 1) {
                    matrix[I*Nx + J][I*Nx + J + 1] = 1;
                    matrix[I*Nx + J][(I + 1)*Nx + J] = 1 + alpha * Hx / lambdas[I][Nx_1];
                } else {
                    matrix[I*Nx + J][I*Nx + J + 1] = top_L / Hx2;
                    matrix[I*Nx + J][(I + 1)*Nx + J] = right_L / Hz2;
                }
            }
        }

        // заполнение диагоналей нижнего блока
        for (int J = 0; J < Nx; J++) {
            matrix[Nz_1*Nx + J][Nz_1*Nx_1 + J] = 1;
            matrix[Nz_1*Nx + J][Nz_1*Nx + J] = 1 + alpha * Hz / lambdas[Nz_1][J];
            matrix[(Nz-2)*Nx + J][Nz_1*Nx + J] = 1 + alpha * Hz / lambdas[Nz_1][J];
        }

        ui->pb->setValue(50);

        prev_Ts = Ts;
        this->resolve_gauss(matrix, Ts, Nx, Nz);
        ui->pb->setValue(70);

        rs.recalc_lambdas(Nx, Nz, Ts, lambdas);
        ui->pb->setValue(90);

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

void MainWindow::resolve_gauss(Doubles2D & matrix, Doubles2D & Ts, int Nx, int Nz) throw (QString) {
    methodGaus gaus(matrix);
    Doubles solution = gaus.calculate();

    for (int I = 0; I < Nz; I++) {
        for (int J = 0; J < Nx; J++) {
            Ts[I][J] = solution[I * Nz + J];
        }
    }
}
