#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "resources.h"

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

    Doubles2D lambdas;
    Resources::init_lambdas(U0, Nx, Nz, lambdas);

    Doubles2D Ts, prev_Ts;
    Resources::init_Ts(Nx, Nz, Ts);
    Resources::init_Ts(Nx, Nz, prev_Ts);

    //TODO: solve task
    ui->error->setText("");
}
