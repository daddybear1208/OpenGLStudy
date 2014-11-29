#include "MainWindow.h"
MainWindow::MainWindow()
{
	setupUi(this);
	pcViewer_ = new PCViewer(this);
	setCentralWidget(pcViewer_);

	connect(actionQuit, SIGNAL(triggered()), this, SLOT(quit()));
}

MainWindow::~MainWindow()
{
	
}

void MainWindow::openFile()
{

}

void MainWindow::closeFile()
{

}

void MainWindow::quit()
{
	this->close();
}

void MainWindow::closeEvent(QCloseEvent* event)
{
	qApp->closeAllWindows();
}

