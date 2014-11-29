#ifndef PCVIEWER_H
#define PCVIEWER_H

#include <QGLViewer/qglviewer.h>

class PCViewer :public QGLViewer
{
public:
	PCViewer(QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags flags = 0);
	PCViewer(QGLContext *context, QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags flags = 0);
	PCViewer(const QGLFormat& format, QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags flags = 0);

	~PCViewer();

protected:
	virtual void draw();
	virtual void init();

private:

};

#endif //! PCIEWER_H