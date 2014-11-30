#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <eigen3/Eigen/Core>

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>

#include "ChunkedArray.h"

#define SQR(x) ((x)*(x))

class PointCloud
{
public:
	PointCloud()
		:point_size_(1.0),
		scale_(1.0),
		shift_(0.0, 0.0, 0.0)
	{
		points_ = new ChunkedArray < 3, float >();
	}

	~PointCloud()
	{
		delete points_;
	}

	// get points count
	inline unsigned getPointsCount() const { return points_->size(); }
	// get point size
	inline const float& getPointSize() const { return point_size_; }
	// set point size
	inline void setPointSize(float size) { point_size_ = size; }
	// add point
	inline void addPoint(const Eigen::Vector3f& pt){ points_->addElement(pt.data()); }
	// get the bounding box
	inline void getBoundingBox(float* minCorner, float* maxCorner) const 
	{
		memcpy(minCorner, points_->minValsPtr(), 3 * sizeof(float));
		memcpy(maxCorner, points_->maxValsPtr(), 3 * sizeof(float));
	}
	// get the real center of point cloud, regardless of the shift and scale
	inline Eigen::Vector3f getCenter() const
	{
		return Eigen::Vector3f((points_->minVal(0) + points_->maxVal(0)) / 2.0,
			(points_->minVal(1) + points_->maxVal(1)) / 2.0,
			(points_->minVal(2) + points_->maxVal(2)) / 2.0);
	}
	// get the radius of the bounding sphere
	inline float getRadius() const 
	{
		return std::sqrt(SQR(points_->maxVal(0)-points_->minVal(0))+
			SQR(points_->maxVal(1) - points_->minVal(1))+
			SQR(points_->maxVal(2) - points_->minVal(2)));
	}
	// draw the point cloud
	void draw()
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		unsigned int chunksCount = points_->chunksCount();
		for (unsigned int i = 0; i < chunksCount; ++i)
		{
			glVertexPointer(3, GL_FLOAT, 3 * sizeof(float), points_->chunkStartPtr(i));
			glDrawArrays(GL_POINTS, 0, (GLsizei)(points_->elementCountOfChunk(i)));
		}
		glDisableClientState(GL_VERTEX_ARRAY);
		
		//glBegin(GL_POINTS);
		//for (unsigned int i = 0; i < points_->size(); ++i)
		//{
		//	glVertex3fv(points_->value(i));
		//}
		//glEnd();

	}

private:
	ChunkedArray<3, float>* points_;
	Eigen::Vector3f shift_;
	float scale_;
	float point_size_;

};


#endif // POINTCLOUD_H
