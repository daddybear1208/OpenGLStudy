#include "MainWindow.h"
#include "AsciiPCReader.h"
#include "LasPCReader.h"


MainWindow::MainWindow()
{
	setupUi(this);
	pcViewer_ = new PCViewer(this);
	setCentralWidget(pcViewer_);

	connect(actionQuit, SIGNAL(triggered()), this, SLOT(quit()));
	connect(this->actionOpen, SIGNAL(triggered()), this, SLOT(openFile()));
	connect(this->actionClose, SIGNAL(triggered()), this, SLOT(closeFile()));
}

MainWindow::~MainWindow()
{
	delete pcViewer_;
}

void MainWindow::openFile()
{
	QString pcFile = QFileDialog::getOpenFileName(this, "Open", "", tr("Point Cloud (*.txt *.las)"));
	if (pcFile.isEmpty()) return;
	pcViewer_->clearPointCloud();
	PointCloud* pointCloud = new PointCloud;
	PCReader::readPCFile(pointCloud, pcFile);
	pcViewer_->setPointCloud(pointCloud);
}

void MainWindow::closeFile()
{
	pcViewer_->clearPointCloud();
}

void MainWindow::quit()
{
	this->close();
}

void MainWindow::closeEvent(QCloseEvent* event)
{
	qApp->closeAllWindows();
}

