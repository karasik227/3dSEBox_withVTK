#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <array>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->GRAPH_2D->hide();

    vtkNew<vtkNamedColors> colors;

    std::array<unsigned char, 4> bkg{{26, 51, 102, 255}};
    colors->SetColor("BkgColor", bkg.data());

    vtkNew<vtkCylinderSource> cylinder;
    cylinder->SetResolution(8);


    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(cylinder->GetOutputPort());

    vtkNew<vtkActor> cylinderActor;
    cylinderActor->SetMapper(mapper);
    cylinderActor->GetProperty()->SetColor(
        colors->GetColor4d("Tomato").GetData());
    cylinderActor->RotateX(30.0);
    cylinderActor->RotateY(-45.0);

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(cylinderActor);
    renderer->SetBackground(1.0, 1.0, 1.0);


    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);

    ui->qvtkWidget_3D_MODEL->setRenderWindow(renderWindow);
    ui->qvtkWidget_3D_MODEL->update();

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_2D_clicked()
{
    ui->qvtkWidget_GRAPH->hide();
    if (graph_2d_exists)
    {
        ui->GRAPH_2D->show();
    }
    else
    {
    series = new QLineSeries();
    QChart *chart = new QChart();


    // Создаем виджет для отображения графика

    //series->setName("рассчет");
    //series->setColor(Qt::blue);
    series->setPen(QPen(Qt::blue, 2, Qt::SolidLine));
    //series->setPointsVisible(true);
    chart->addSeries(series);
    chart->createDefaultAxes();
    //chart->setTitle("График данных");
    //chart->setTitleFont(QFont("Times New Roman", 14, QFont::Bold));
    chart->setAnimationOptions(QChart::AllAnimations);

    // Настройка осей
    QValueAxis *axisX = new QValueAxis;
    axisX->setTitleText("Частота, Гц");
    axisX->setTitleFont(QFont("Times New Roman", 12, QFont::Bold));
    axisX->setLabelsFont(QFont("Arial", 10));
    axisX->setLabelFormat("%d");
    axisX->setTickCount(10);
    chart->setAxisX(axisX, series);

    QValueAxis *axisY = new QValueAxis;
    axisY->setTitleText("ЭЭ, дБ");
    axisY->setTitleFont(QFont("Times New Roman", 12, QFont::Bold));
    axisY->setLabelsFont(QFont("Arial", 10));
    axisY->setLabelFormat("%.2f");
    axisY->setTickCount(10);
    chart->setAxisY(axisY, series);


    ui->GRAPH_2D->setChart(chart);

    ui->GRAPH_2D->setRenderHint(QPainter::Antialiasing);
    chart->legend()->hide();
    ui->GRAPH_2D->show();

    graph_2d_exists = true;
    }
}


void MainWindow::on_pushButton_3D_clicked()
{
    ui->qvtkWidget_GRAPH->show();
    ui->GRAPH_2D->hide();
}

