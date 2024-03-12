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

void MainWindow::minimizeCostFunction(int firstSample, int batchSize){
    QDebug debug = qDebug();

    QList<QList<QList<double>>> wGradient;
    QList<QList<double>> bGradient;
    double Co;
    for (int L = 0; L < layersNum; ++L) {
        wGradient.append(QList<QList<double>>());
        bGradient.append(QList<double>());
        for (int j = 0; j < neuronsNum[L+1]; ++j) {
            wGradient[L].append(QList<double>());
            for(int k = 0; k < neuronsNum[L]; ++k){
                wGradient[L][j].append(0);
            }
            bGradient[L].append(0);
        }
    }

    //Compute gradient descent's vector--------------------------------------
    for (int image = firstSample; image < firstSample + batchSize; ++image) {
//        debug << "Sample#" << image << Qt::endl
//              << Qt::endl;

        //Set objective outputs
        QList<int> y;//Objective outputs
        for (int i = 0; i < neuronsNum[layersNum]; ++i) {
            if(i == (int)labelsFile.at(image)){
                y.append(1);
            }else{
                y.append(0);
            }
        }
//        debug << "Objective :" << y << Qt::endl
//              << Qt::endl;

        //Set inputs
        for (int i = 0; i < neuronsNum[0]; ++i) {
            a[0][i] = (double)imagesFile.at(image).at(i) / 256;
        }

        //Compute outputs
        QList<QList<double>> z;//weighted sum
        z.append(QList<double>());
        for (int endLayer = 1; endLayer < layersNum + 1; ++endLayer) {
            z.append(QList<double>());
            for (int startLayer = 0; startLayer < neuronsNum[endLayer]; ++startLayer) {
                double n = 0;
                for (int k = 0; k < neuronsNum[endLayer-1]; ++k) {
                    n += w[endLayer-1][startLayer][k] * a[endLayer-1][k];
                }
                n += b[endLayer-1][startLayer];
                z[endLayer].append(n);
                a[endLayer][startLayer] = (tanh(n) + 1) / 2;
            }
//            debug << "Layer " << endLayer+1 << ":"
//                  << a[endLayer] << Qt::endl
//                  << neuronsNum[endLayer] << "x" << neuronsNum[endLayer-1] << Qt::endl;
        }
//        debug << Qt::endl;

        //Compute cost function
        Co = 0;//cost function
        QList<double> cost;
        for (int i = 0; i < neuronsNum[layersNum]; ++i) {
            cost.append(pow((a[layersNum][i] - y[i]), 2));
            Co += cost[i];
        }
        Co /= neuronsNum[layersNum];
//        debug << "Costs :" << cost << Qt::endl
//              << "Cost function: " << Co[image] << Qt::endl
//              << Qt::endl;

        //Compute gradient descent's derivatives
        QList<QList<double>> cost_z;//da/dz * dCo/da
        for (int var = 0; var < layersNum; ++var) {
            cost_z.append(QList<double>());
        }
        int L = layersNum;
        for (int k = 0; k < neuronsNum[L]; ++k) {
            double cost_a = 2 * (a[L][k] - y[k]);
            double a_z = pow(sech(z[L][k]), 2);
            cost_z[L - 1].append(a_z * cost_a);

            bGradient[L - 1][k] += cost_z[L - 1][k];
            for (int j = 0; j < neuronsNum[L - 1]; ++j) {
                wGradient[L - 1][k][j] += a[L][k] * cost_z[L - 1][k];
            }
        }
        for (int L = layersNum - 1; L > 0; --L) {
            for (int k = 0; k < neuronsNum[L]; ++k) {
                double cost_a = 0;
                for (int j = 0; j < neuronsNum[L + 1]; ++j) {
                    double z_a = w[L][j][k];
                    cost_a += z_a * cost_z[L][j];
                }
                double a_z = pow(sech(z[L][k]), 2);
                cost_z[L - 1].append(a_z * cost_a);

                bGradient[L - 1][k] += cost_z[L - 1][k];
                for (int j = 0; j < neuronsNum[L - 1]; ++j) {
                    wGradient[L - 1][k][j] += a[L][k] * cost_z[L - 1][k];
                }
            }
        }
//        debug << "Gradient of w :" << Qt::endl
//              << wGradient << Qt::endl
//              << "Gradient of b :" << Qt::endl
//              << bGradient << Qt::endl
//              << Qt::endl;
    }

    for (int L = 0; L < layersNum; ++L) {
        for (int j = 0; j < neuronsNum[L+1]; ++j) {
            for(int k = 0; k < neuronsNum[L]; ++k){
                wGradient[L][j][k] /= batchSize;
                w[L][j][k] -= 0.01 * wGradient[L][j][k];
            }
            bGradient[L][j] /= batchSize;
            b[L][j] -= 0.01 * bGradient[L][j];
        }
    }
//    debug << "Average :" << Qt::endl
//          << "Gradient of w :" << Qt::endl
//          << wGradient << Qt::endl
//          << "Gradient of b :" << Qt::endl
//          << bGradient << Qt::endl
//          << Qt::endl;
    debug << "Cost#" << firstSample / batchSize << Co;
}

double MainWindow::sech(double x) {
    return 1.0 / cosh(x);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QDebug debug = qDebug();

    for (int L = 0; L < layersNum + 1; ++L) {
        a.append(QList<double>());
        for (int j = 0; j < neuronsNum[L]; ++j) {
            a[L].append(0);
        }
    }

    for (int L = 0; L < layersNum; ++L) {
        w.append(QList<QList<double>>());
        b.append(QList<double>());
        for (int j = 0; j < neuronsNum[L+1]; ++j) {
            w[L].append(QList<double>());
//            debug << "[" << L << "][" << j << "]" << Qt::endl;
//            debug << "Weights : ";
            for(int k = 0; k < neuronsNum[L]; ++k){
                w[L][j].append(randomGenerator.generateDouble() * 2 - 1);
//                debug << w[L][j][k] << " ";
            }
            b[L].append(randomGenerator.generateDouble() / 2);
//            debug << "Bias : " << b[L][j] << Qt::endl;
        }
//        debug << Qt::endl;
    }


    QString filename = QString::fromStdString("C:/Users/david/Documents/QtCreator/Projects/NeuralNetwork/build-NeuralNetwork-Desktop_Qt_5_15_2_MinGW_32_bit-Debug/train-images-idx3-ubyte/train-images-idx3-ubyte");
    QString label_filename = QString::fromStdString("C:/Users/david/Documents/QtCreator/Projects/NeuralNetwork/build-NeuralNetwork-Desktop_Qt_5_15_2_MinGW_32_bit-Debug/train-labels-idx1-ubyte/train-labels-idx1-ubyte");

    imagesFile = readIDX3UByteFile(filename);
    labelsFile = readLabelFile(label_filename);

//    int random = QRandomGenerator::global()->bounded(60001);
//    int exemple = 0;
//    debug << "Random exemple : " << exemple << Qt::endl;
//    for (int i = 0; i < (int)imagesFile.at(exemple).size(); ++i) {
//        int pixel = (int)imagesFile.at(exemple).at(i);
//        debug << pixel;
//        int digits = 1;
//        while(pixel >= 10){
//            pixel /= 10;
//            ++digits;
//        }
//        for(int j = 0; j < 3 - digits; ++j){
//            debug << "";
//        }
//        if((i + 1) % 28 == 0){
//            debug << Qt::endl;
//        }
//    }
//    if (!labelsFile.empty()) {
//        debug << "Label : " << labelsFile.at(exemple);
//    } else {
//        debug << "Labels vector is empty or does not contain any elements.";
//    }
//    debug << Qt::endl;

    for (int i = 0; i < 600; ++i) {
        minimizeCostFunction(100 * i, 100);
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

