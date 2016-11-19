#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "resources.h"

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
    void resolve_gauss(Doubles2D & matrix, Doubles2D & Ts);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
