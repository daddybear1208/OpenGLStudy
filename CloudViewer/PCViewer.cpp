#include "PCViewer.h"

PCViewer::PCViewer(QWidget* parent /*= 0*/, const QGLWidget* shareWidget /*= 0*/, Qt::WindowFlags flags /*= 0*/)
	:QGLViewer(parent, shareWidget, flags),
	point_cloud_(0)
{
	
}

PCViewer::PCViewer(QGLContext *context, QWidget* parent /*= 0*/, const QGLWidget* shareWidget /*= 0*/, Qt::WindowFlags flags /*= 0*/)
	: QGLViewer(context, parent, shareWidget, flags),
	point_cloud_(0)
{

}

PCViewer::PCViewer(const QGLFormat& format, QWidget* parent /*= 0*/, const QGLWidget* shareWidget /*= 0*/, Qt::WindowFlags flags /*= 0*/)
	: QGLViewer(format, parent, shareWidget, flags),
	point_cloud_(0)
{

}

PCViewer::~PCViewer()
{
	if (point_cloud_)
	{
		delete point_cloud_;
	}
}

void PCViewer::draw()
{
	if (point_cloud_)
		point_cloud_->draw();
}

void PCViewer::init()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glDisable(GL_LIGHTING);
}

void PCViewer::setPointCloud(PointCloud* pc)
{
	point_cloud_ = pc;
	Eigen::Vector3f center = pc->getCenter();
	setSceneCenter(qglviewer::Vec(center[0], center[1], center[2]));
	setSceneRadius(pc->getRadius());
	showEntireScene();
}

void PCViewer::clearPointCloud()
{
	if (point_cloud_)
	{
		delete point_cloud_;
	}
	point_cloud_ = 0;
	setSceneCenter(qglviewer::Vec(0.0, 0.0, 0.0));
	setSceneRadius(100);
}



