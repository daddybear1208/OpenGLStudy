#ifndef MAINWINDOW_H
#define	MAINWINDOW_H

#include <QtGui>
#include <QMainWindow>
#include <QFileDialog>

#include "PCViewer.h"
#include "ui_MainWindow.h"

class MainWindow :public QMainWindow, public Ui::MainWindow
{
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();

protected:
	void closeEvent(QCloseEvent* event);

protected slots:
	void openFile();
	void closeFile();
	void quit();



private:
	PCViewer* pcViewer_;

};

#endif // !MAINWINDOW_H



