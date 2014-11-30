#ifndef PCVIEWER_H
#define PCVIEWER_H

#include <QGLViewer/qglviewer.h>

#include "PointCloud.h"

class PCViewer :public QGLViewer
{
	//Q_OBJECT
public:
	PCViewer(QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags flags = 0);
	PCViewer(QGLContext *context, QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags flags = 0);
	PCViewer(const QGLFormat& format, QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags flags = 0);

	~PCViewer();

public slots:
	void setPointCloud(PointCloud* pc);
	void clearPointCloud();

protected:
	virtual void draw();
	virtual void init();

private:
	PointCloud* point_cloud_;
};

#endif //! PCIEWER_H