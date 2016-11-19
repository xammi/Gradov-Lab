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
    try {
        Doubles2D Ts;
        this->recalculate_action(Ts);
    }
    catch (QString error) {
        ui->error->setText(error);
    }
}

void MainWindow::recalculate_action(Doubles2D & Ts) throw (QString) {
    double A = ui->A->value(), B = ui->B->value();
    double Ft = ui->Ft->value(), F0 = ui->F0->value();
    double U0 = ui->U0->value(), alpha = ui->alpha->value();
    int Nx = ui->Nx->value(), Nz = ui->Nz->value();

    if (Nx == 0) {
        ui->error->setText("Количество узлов сетки по X равно 0");
        return;
    }
    if (Nz == 0) {
        ui->error->setText("Количество узлов сетки по Z равно 0");
        return;
    }

    double Hx = A / Nx, Hz = B / Nz;
    double Hx2 = Hx * Hx, Hz2 = Hz * Hz;
    int dim = Nx * Nz;

    Doubles2D lambdas;
    Resources::init_lambdas(U0, Nx, Nz, lambdas);

    Doubles2D prev_Ts;
    Resources::init_Ts(Nx, Nz, Ts);

    Doubles2D matrix;
    Resources::init_matrix(Nx, Nz, matrix);

    while (Resources::if_stop_iterations(Nx, Nz, Ts, prev_Ts)) {

        // вычисление 1-го краевого условия (T[0][J] - T[1][J] = K1[J])
        for (int J = 0; J < Nx; J++) {
            double K1 = Hz * Ft / lambdas[0][J];
            matrix[dim][J] = K1;
        }

        // вычисление 2-го краевого условия (K1[I] * T[I][0] - T[I][1] = K2[I])
        for (int I = 0; I < Nz; I++) {
            double K1 = 1 - alpha * Hx / lambdas[I][0];
            double K2 = -alpha * Hx / lambdas[I][0] * U0;
            matrix[dim][Nx + Nx * I] = K2;
        }

        // вычисление 3-го краевого условия (T[I][Nx - 1] - K1[I] * T[I][Nx] = K2[I])
        for (int I = 0; I < Nz; I++) {
            double K1 = 1 + alpha * Hx / lambdas[I][Nx];
            double K2 = -alpha * Hx / lambdas[I][Nx] * U0;
            matrix[dim][2*Nx - 1 + Nx * I] = K2;
        }

        // вычисление 4-го краевого условия (T[Nz - 1][J] - K1[J] * T[Nz][J] = K2[J])
        for (int J = 0; J < Nx; J++) {
            double K1 = 1 + alpha * Hz / lambdas[Nz][J];
            double K2 = -alpha * Hz / lambdas[Nz][J] * U0;
            matrix[dim][dim - Nx + J] = K2;
        }

        // заполнение правой части
        for (int I = 1, J = 1; I < Nz - 1; I++, J++) {
            matrix[dim][I * Nx + J] = Resources::calc_f(F0, Hz * I, Hx * J);
        }

        // заполнение диагоналей
        for (int I = 1, J = 1; I < Nz - 1; I++, J++) {
            double bottom_L = (lambdas[I][J] + lambdas[I][J - 1]) / 2;
            double top_L = (lambdas[I][J] + lambdas[I][J + 1]) / 2;
            double left_L = (lambdas[I][J] + lambdas[I - 1][J]) / 2;
            double right_L = (lambdas[I][J] + lambdas[I + 1][J]) / 2;

            matrix[I][J - 1] = left_L / Hx2;
            matrix[I][J + 1] = top_L / Hx2;

            matrix[I][J] = (top_L + bottom_L) / Hx2 + (left_L + right_L) / Hz2;

            matrix[I - 1][J] = left_L / Hz2;
            matrix[I + 1][J] = right_L / Hz2;
        }

        prev_Ts = Ts;
        this->resolve_gauss(matrix, Ts);
        Resources::recalc_lambdas(Nx, Nz, Ts, lambdas);
    }
}

void resolve_gauss(Doubles2D & matrix, Doubles2D & Ts) {

}
