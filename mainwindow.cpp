#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->btn_recalc, SIGNAL(pressed()), SLOT(recalculate()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::recalculate() {
    double A = ui->A->value(), B = ui->B->value();
    double Nx = ui->Nx->value(), Nz = ui->Nz->value();
    double Ft = ui->Ft->value(), F0 = ui->F0->value();
    double U0 = ui->U0->value(), alpha = ui->alpha->value();

    //TODO: solve task
}
