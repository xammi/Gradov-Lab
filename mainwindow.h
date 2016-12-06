#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "resources.h"
#include "methodgaus.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public slots:
    void recalculate();

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void recalculate_action(Doubles2D & Ts) throw (QString);
    void resolve_gauss(Doubles2D & matrix, Doubles2D & Ts, int Nx, int Nz) throw (QString) ;
    void view_result(Doubles2D & Ts);

private:
    Ui::MainWindow *ui;
    Resources rs;
    methodGaus* gaus;
};

#endif // MAINWINDOW_H
