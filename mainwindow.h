#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include <QRandomGenerator>
#include <QFile>
#include <QDataStream>
#include <QtCharts/QtCharts>

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

    QList<int> neuronsPerLayerList = {728, 16, 16, 10};//neurons number per layer
    QList<QList<double>> a;//neurons values
    QList<QList<QList<double>>> w;//weights
    QList<QList<double>> b;//biases

    int batchSize = 100;
    double r = 0.1;//Learning rate

    QList<QList<unsigned char>> imagesFile;
    QList<unsigned char> labelsFile;
    
    void minimizeCostFunction(int firstSample, int batchSize, double &costAverage, int &accuracy);
    double sech(double x);
};
#endif // MAINWINDOW_H
