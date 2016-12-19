#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "resources.h"
#include "gauss.h"

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

    void straight_method(Doubles2D & Ts) throw (QString);
    void relaxation_method(Doubles2D & Ts) throw (QString);
    void view_result(Doubles2D & Ts);

private:
    Ui::MainWindow *ui;
    Resources rs;
    GaussResolver* gauss;
};

#endif // MAINWINDOW_H
