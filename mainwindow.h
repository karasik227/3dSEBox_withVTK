#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_2D_clicked();

    void on_pushButton_3D_clicked();

private:
    Ui::MainWindow *ui;
    QLineSeries *series;
    QChartView *chartView;
    bool graph_2d_exists = false;
};
#endif // MAINWINDOW_H
