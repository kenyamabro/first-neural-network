#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <fstream>
#include <vector>

QList<QList<unsigned char>> readIDX3UByteFile(const QString& filename) {
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "Failed to open the IDX3-UBYTE file.";
        return {};
    }

    QDataStream in(&file);
    in.setByteOrder(QDataStream::BigEndian);

    // Read the IDX3-UBYTE file header
    char magicNumber[4];
    in.readRawData(magicNumber, 4);

    char numImagesBytes[4];
    char numRowsBytes[4];
    char numColsBytes[4];

    in.readRawData(numImagesBytes, 4);
    in.readRawData(numRowsBytes, 4);
    in.readRawData(numColsBytes, 4);

    // Convert the header information from big-endian to native endianness
    int numImages = (static_cast<unsigned char>(numImagesBytes[0]) << 24) |
                    (static_cast<unsigned char>(numImagesBytes[1]) << 16) |
                    (static_cast<unsigned char>(numImagesBytes[2]) << 8) |
                    static_cast<unsigned char>(numImagesBytes[3]);
    int numRows = (static_cast<unsigned char>(numRowsBytes[0]) << 24) |
                  (static_cast<unsigned char>(numRowsBytes[1]) << 16) |
                  (static_cast<unsigned char>(numRowsBytes[2]) << 8) |
                  static_cast<unsigned char>(numRowsBytes[3]);
    int numCols = (static_cast<unsigned char>(numColsBytes[0]) << 24) |
                  (static_cast<unsigned char>(numColsBytes[1]) << 16) |
                  (static_cast<unsigned char>(numColsBytes[2]) << 8) |
                  static_cast<unsigned char>(numColsBytes[3]);

    // Initialize a QList to store the images
    QList<QList<unsigned char>> images;

    for (int i = 0; i < numImages; ++i) {
        // Read each image as a QList of bytes
        QList<unsigned char> image;
        for (int j = 0; j < numRows * numCols; ++j) {
            unsigned char pixel;
            in >> pixel;
            image.append(pixel);
        }
        images.append(image);
    }

    file.close();

    return images;
}

QList<unsigned char> readLabelFile(const QString& filename) {
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "Failed to open the IDX3-UBYTE file.";
        return {};
    }

    QDataStream in(&file);
    in.setByteOrder(QDataStream::BigEndian);

    // Read the IDX3-UBYTE file header
    char magicNumber[4];
    char numImagesBytes[4];

    in.readRawData(magicNumber, 4);
    in.readRawData(numImagesBytes, 4);

    // Convert the header information from big-endian to native endianness
    int numImages = (static_cast<unsigned char>(numImagesBytes[0]) << 24) |
                    (static_cast<unsigned char>(numImagesBytes[1]) << 16) |
                    (static_cast<unsigned char>(numImagesBytes[2]) << 8) |
                    static_cast<unsigned char>(numImagesBytes[3]);

    // Initialize a QList to store the labels
    QList<unsigned char> labels;

    for (int i = 0; i < numImages; ++i) {
        // Read each label as a single byte
        unsigned char label;
        in >> label;
        labels.append(label);
    }

    file.close();

    return labels;
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QDebug debug = qDebug();

    this->setCentralWidget(new QWidget());
    auto layout = new QHBoxLayout();
    this->centralWidget()->setLayout(layout);

    auto initializerLayout = new QFormLayout();
    layout->addLayout(initializerLayout, 1);

    auto initializeButton = new QPushButton("Initialize");
    initializerLayout->addRow(initializeButton);

    int batchSize = 100;

    auto costSeries = new QLineSeries;
    auto acuracySeries = new QLineSeries;

    costSeries->setName("Cost");
    acuracySeries->setName("Acuracy");

    QChart* chart = new QChart();
    chart->setTitle("Cost and acuracy");
    QFont font;
    font.setPointSize(12);
    font.setBold(true);
    chart->setTitleFont(font);

    costSeries->setPointsVisible(true);
    acuracySeries->setPointsVisible(true);
    chart->addSeries(costSeries);
    chart->addSeries(acuracySeries);

    chart->createDefaultAxes();

    chart->axes(Qt::Horizontal).constFirst()->setTitleText("Trainings");
    chart->axes(Qt::Vertical).constFirst()->setTitleText("Decimal");

    auto xAxis = qobject_cast<QValueAxis*>(chart->axes(Qt::Horizontal).at(0));
    auto yAxis = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).at(0));

    xAxis->setTickType(QValueAxis::TicksFixed);
    yAxis->setTickType(QValueAxis::TicksFixed);

    xAxis->setLabelFormat("%g");
    yAxis->setLabelFormat("%g");

    yAxis->setMin(0);
    yAxis->setMax(1);

    auto chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    layout->addWidget(chartView, 5);

    QObject::connect(initializeButton, &QPushButton::clicked, this, [=]{
        for (int L = 0; L < neuronsNumPerLayer.size(); ++L) {
            a.append(QList<double>());
            for (int j = 0; j < neuronsNumPerLayer[L]; ++j) {
                a[L].append(0);
            }
        }

        for (int L = 0; L < neuronsNumPerLayer.size() - 1; ++L) {
            w.append(QList<QList<double>>());
            b.append(QList<double>());
            for (int j = 0; j < neuronsNumPerLayer[L+1]; ++j) {
                w[L].append(QList<double>());
                for(int k = 0; k < neuronsNumPerLayer[L]; ++k){
                    w[L][j].append(QRandomGenerator::global()->generateDouble() * 2 - 1);
                }
                b[L].append(QRandomGenerator::global()->generateDouble() / 2);
            }
        }

        QString filename = QString::fromStdString("C:/Users/david/Documents/QtCreator/Projects/NeuralNetwork/build-NeuralNetwork-Desktop_Qt_5_15_2_MinGW_32_bit-Debug/train-images-idx3-ubyte/train-images-idx3-ubyte");
        QString label_filename = QString::fromStdString("C:/Users/david/Documents/QtCreator/Projects/NeuralNetwork/build-NeuralNetwork-Desktop_Qt_5_15_2_MinGW_32_bit-Debug/train-labels-idx1-ubyte/train-labels-idx1-ubyte");

        imagesFile = readIDX3UByteFile(filename);
        labelsFile = readLabelFile(label_filename);

        for (int x = 0; x < 600; ++x) {
            double cost = 0;
            int acuracy = 0;

            minimizeCostFunction(batchSize * x, batchSize, cost, acuracy);
            emit batchTrained(x, cost, acuracy);

            QCoreApplication::processEvents();
        }
    });

    connect(this, &MainWindow::batchTrained, this, [=](int x, double cost, int acuracy) {
        costSeries->append(QPointF(x, cost));
        acuracySeries->append(QPointF(x, (double)acuracy / batchSize));

        xAxis->setMax(x);
    });
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::minimizeCostFunction(int firstSample, int batchSize, double &costAverage, int &accuracy){
    QDebug debug = qDebug();

    QList<QList<QList<double>>> wGradient;
    QList<QList<double>> bGradient;
    for (int L = 0; L < neuronsNumPerLayer.size() - 1; ++L) {
        wGradient.append(QList<QList<double>>());
        bGradient.append(QList<double>());
        for (int j = 0; j < neuronsNumPerLayer[L+1]; ++j) {
            wGradient[L].append(QList<double>());
            for(int k = 0; k < neuronsNumPerLayer[L]; ++k){
                wGradient[L][j].append(0);
            }
            bGradient[L].append(0);
        }
    }

    //Compute gradient descent's vector--------------------------------------
    for (int image = firstSample; image < firstSample + batchSize; ++image) {
        //Set objective outputs (y)
        QList<int> y;//Objective outputs
        for (int i = 0; i < neuronsNumPerLayer[neuronsNumPerLayer.size() - 1]; ++i) {
            if(i == (int)labelsFile.at(image)){
                y.append(1);
            }else{
                y.append(0);
            }
        }

        //Set inputs (a[0])
        for (int i = 0; i < neuronsNumPerLayer[0]; ++i) {
            a[0][i] = (double)imagesFile.at(image).at(i) / 256;
        }

        //Compute outputs (a)
        QList<QList<double>> z;//weighted sum
        z.append(QList<double>());
        for (int endLayer = 1; endLayer < neuronsNumPerLayer.size(); ++endLayer) {
            z.append(QList<double>());
            for (int startLayer = 0; startLayer < neuronsNumPerLayer[endLayer]; ++startLayer) {
                double n = 0;
                for (int k = 0; k < neuronsNumPerLayer[endLayer-1]; ++k) {
                    n += w[endLayer-1][startLayer][k] * a[endLayer-1][k];
                }
                n += b[endLayer-1][startLayer];
                z[endLayer].append(n);
                a[endLayer][startLayer] = (tanh(n) + 1) / 2;
            }
        }

        //Compute cost function
        double Co = 0;//cost function
        QList<double> cost;
        double biggestCost = 0;
        int biggestIndex = -1;
        for (int i = 0; i < neuronsNumPerLayer[neuronsNumPerLayer.size() - 1]; ++i) {
            cost.append(pow((a[neuronsNumPerLayer.size() - 1][i] - y[i]), 2));
            Co += cost[i];
            if(cost[i] > biggestCost){
                biggestCost = cost[i];
                biggestIndex = i;
            }
        }
        if(biggestIndex == (int)labelsFile.at(image)) ++accuracy;
        Co /= neuronsNumPerLayer[neuronsNumPerLayer.size() - 1];
        costAverage += Co;

        //Compute gradient descent's derivatives (back propagation)
        QList<QList<double>> cost_z;//da/dz * dCo/da
        for (int var = 0; var < neuronsNumPerLayer.size() - 1; ++var) {
            cost_z.append(QList<double>());
        }
        int L = neuronsNumPerLayer.size() - 1;
        for (int k = 0; k < neuronsNumPerLayer[L]; ++k) {
            double cost_a = 2 * (a[L][k] - y[k]);//dCo/da
            double a_z = pow(sech(z[L][k]), 2);//da/dz
            cost_z[L - 1].append(a_z * cost_a);

            bGradient[L - 1][k] += cost_z[L - 1][k];//dCo/db
            for (int j = 0; j < neuronsNumPerLayer[L - 1]; ++j) {
                wGradient[L - 1][k][j] += a[L][k] * cost_z[L - 1][k];//dCo/dw
            }
        }
        for (int L = neuronsNumPerLayer.size() - 2; L > 0; --L) {
            for (int k = 0; k < neuronsNumPerLayer[L]; ++k) {
                double cost_a = 0;//dCo/da
                for (int j = 0; j < neuronsNumPerLayer[L + 1]; ++j) {
                    double z_a = w[L][j][k];//dz(L)/da(L-1)
                    cost_a += z_a * cost_z[L][j];
                }
                double a_z = pow(sech(z[L][k]), 2);//da/dz
                cost_z[L - 1].append(a_z * cost_a);

                bGradient[L - 1][k] += cost_z[L - 1][k];//dCo/db
                for (int j = 0; j < neuronsNumPerLayer[L - 1]; ++j) {
                    wGradient[L - 1][k][j] += a[L][k] * cost_z[L - 1][k];//dCo/dw
                }
            }
        }
    }
    costAverage /= batchSize;

    for (int L = 0; L < neuronsNumPerLayer.size() - 1; ++L) {
        for (int j = 0; j < neuronsNumPerLayer[L+1]; ++j) {
            for(int k = 0; k < neuronsNumPerLayer[L]; ++k){
                wGradient[L][j][k] /= batchSize;
                w[L][j][k] -= 0.01 * wGradient[L][j][k];
            }
            bGradient[L][j] /= batchSize;
            b[L][j] -= 0.01 * bGradient[L][j];
        }
    }
    debug << firstSample / batchSize << "# Cost:" << costAverage
          << "Accuracy:" << accuracy << "/" << batchSize << " ";
}

double MainWindow::sech(double x) {
    return 1.0 / cosh(x);
}
