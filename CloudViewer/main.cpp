#include <iostream>
#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>
#include <QtGui>
#include <QApplication>

#include "MainWindow.h"


int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	MainWindow mainWindow;
	mainWindow.show();
	return app.exec();
}

