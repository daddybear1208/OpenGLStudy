#include "PCViewer.h"

PCViewer::PCViewer(QWidget* parent /*= 0*/, const QGLWidget* shareWidget /*= 0*/, Qt::WindowFlags flags /*= 0*/)
	:QGLViewer(parent, shareWidget, flags)
{

}

PCViewer::PCViewer(QGLContext *context, QWidget* parent /*= 0*/, const QGLWidget* shareWidget /*= 0*/, Qt::WindowFlags flags /*= 0*/)
	: QGLViewer(context, parent, shareWidget, flags)
{

}

PCViewer::PCViewer(const QGLFormat& format, QWidget* parent /*= 0*/, const QGLWidget* shareWidget /*= 0*/, Qt::WindowFlags flags /*= 0*/)
	: QGLViewer(format, parent, shareWidget, flags)
{

}

PCViewer::~PCViewer()
{

}

void PCViewer::draw()
{

}

void PCViewer::init()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
}



