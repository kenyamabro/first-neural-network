#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include <QRandomGenerator>
#include <QFile>
#include <QDataStream>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

    QList<int> neuronsNum = {728, 16, 16, 10};//neurons number per layer
    int layersNum = neuronsNum.size() - 1;//layers index
    QList<QList<double>> a;//neurons values
    QList<QList<QList<double>>> w;//weights
    QList<QList<double>> b;//biases

    QList<QList<unsigned char>> imagesFile;
    QList<unsigned char> labelsFile;
    
    void minimizeCostFunction(int firstSample, int batchSize);

    QRandomGenerator randomGenerator;

    double sech(double x);
};
#endif // MAINWINDOW_H
