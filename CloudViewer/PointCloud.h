#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include "ChunkedArray.h"


#include <eigen3/Eigen/Core>

class PointCloud
{
	PointCloud();
	~PointCloud();

	// get points count
	inline unsigned getPointsCount() const { return points_->size(); }
	// get point size
	inline const float& getPointSize() const { return point_size_; }
	// non-const version of getPointSize
	inline float& getPointSize() { return const_cast<float&>(static_cast<const PointCloud*>(this)->getPointSize()); };
	// allocate a piece of memory for storing the points
	inline bool reserve(unsigned int pointCount){ return points_->reserve(pointCount); }
	// add point
	inline void addPoint(const Eigen::Vector3f& pt){ points_->addElement(pt.data()); }
	// draw point cloud
	void draw();

private:
	ChunkedArray<3, float>* points_;
	Eigen::Vector3f shift_;
	float scale_;
	float point_size_;

};


#endif // POINTCLOUD_H
