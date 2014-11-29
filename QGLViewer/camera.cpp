#include "domUtils.h"
#include "camera.h"
#include "qglviewer.h"
#include "manipulatedCameraFrame.h"

using namespace std;
using namespace qglviewer;

Camera::Camera()
	: frame_(NULL), fieldOfView_(M_PI/4.0f), modelViewMatrixIsUpToDate_(false), projectionMatrixIsUpToDate_(false)
{
	// #CONNECTION# Camera copy constructor
	interpolationKfi_ = new KeyFrameInterpolator;
	// Requires the interpolationKfi_
	setFrame(new ManipulatedCameraFrame());

	// #CONNECTION# All these default values identical in initFromDOMElement.

	// Requires fieldOfView() to define focusDistance()
	setSceneRadius(100.0);

	// Initial value (only scaled after this)
	orthoCoef_ = tan(fieldOfView()/2.0);

	// Also defines the revolveAroundPoint(), which changes orthoCoef_. Requires a frame().
	setSceneCenter(Vec(0.0, 0.0, 0.0));

	// Requires fieldOfView() when called with ORTHOGRAPHIC. Attention to projectionMatrix_ below.
	setType(PERSPECTIVE);

	// #CONNECTION# initFromDOMElement default values
	setZNearCoefficient(0.005f);
	setZClippingCoefficient(sqrt(3.0));

	// Dummy values
	setScreenWidthAndHeight(600, 400);

	// Stereo parameters
	setIODistance(0.062f);
	setPhysicalScreenWidth(0.5f);
	// focusDistance is set from setFieldOfView()

	// #CONNECTION# Camera copy constructor
	for (unsigned short j=0; j<16; ++j)
	{
		modelViewMatrix_[j] = ((j%5 == 0) ? 1.0 : 0.0);
		// #CONNECTION# computeProjectionMatrix() is lazy and assumes 0.0 almost everywhere.
		projectionMatrix_[j] = 0.0;
	}
	computeProjectionMatrix();
}


Camera::~Camera()
{
	delete frame_;
	delete interpolationKfi_;
}


/*! Copy constructor. Performs a deep copy using operator=(). */
Camera::Camera(const Camera& camera)
	: QObject()
{
	// #CONNECTION# Camera constructor
	interpolationKfi_ = new KeyFrameInterpolator;
	// Requires the interpolationKfi_
	setFrame(new ManipulatedCameraFrame());

	for (unsigned short j=0; j<16; ++j)
	{
		modelViewMatrix_[j] = ((j%5 == 0) ? 1.0 : 0.0);
		// #CONNECTION# computeProjectionMatrix() is lazy and assumes 0.0 almost everywhere.
		projectionMatrix_[j] = 0.0;
	}

	(*this)=camera;
}

/*! Equal operator.

 All the parameters of \p camera are copied. The frame() pointer is not modified, but its
 Frame::position() and Frame::orientation() are set to those of \p camera.

 \attention The Camera screenWidth() and screenHeight() are set to those of \p camera. If your
 Camera is associated with a QGLViewer, you should update these value after the call to this method:
 \code
 *(camera()) = otherCamera;
 camera()->setScreenWidthAndHeight(width(), height());
 \endcode
 The same applies to sceneCenter() and sceneRadius(), if needed. */
Camera& Camera::operator=(const Camera& camera)
{
	setScreenWidthAndHeight(camera.screenWidth(), camera.screenHeight());
	setFieldOfView(camera.fieldOfView());
	setSceneRadius(camera.sceneRadius());
	setSceneCenter(camera.sceneCenter());
	setZNearCoefficient(camera.zNearCoefficient());
	setZClippingCoefficient(camera.zClippingCoefficient());
	setType(camera.type());

	// Stereo parameters
	setIODistance(camera.IODistance());
	setFocusDistance(camera.focusDistance());
	setPhysicalScreenWidth(camera.physicalScreenWidth());

	orthoCoef_ = camera.orthoCoef_;
	projectionMatrixIsUpToDate_ = false;

	// frame_ and interpolationKfi_ pointers are not shared.
	frame_->setReferenceFrame(NULL);
	frame_->setPosition(camera.position());
	frame_->setOrientation(camera.orientation());

	interpolationKfi_->resetInterpolation();

	kfi_ = camera.kfi_;

	computeProjectionMatrix();
	computeModelViewMatrix();

	return *this;
}

/*! Sets Camera screenWidth() and screenHeight() (expressed in pixels).

You should not call this method when the Camera is associated with a QGLViewer, since the
latter automatically updates these values when it is resized (hence overwritting your values).

Non-positive dimension are silently replaced by a 1 pixel value to ensure frustrum coherence.

If your Camera is used without a QGLViewer (offscreen rendering, shadow maps), use setAspectRatio()
instead to define the projection matrix. */
void Camera::setScreenWidthAndHeight(int width, int height)
{
	// Prevent negative and zero dimensions that would cause divisions by zero.
	screenWidth_  = width > 0 ? width : 1;
	screenHeight_ = height > 0 ? height : 1;
	projectionMatrixIsUpToDate_ = false;
}

/*! Returns the near clipping plane distance used by the Camera projection matrix.

 The clipping planes' positions depend on the sceneRadius() and sceneCenter() rather than being fixed
 small-enough and large-enough values. A good scene dimension approximation will hence result in an
 optimal precision of the z-buffer.

 The near clipping plane is positioned at a distance equal to zClippingCoefficient() * sceneRadius()
 in front of the sceneCenter():
 \code
 zNear = distanceToSceneCenter() - zClippingCoefficient()*sceneRadius();
 \endcode

 In order to prevent negative or too small zNear() values (which would degrade the z precision),
 zNearCoefficient() is used when the Camera is inside the sceneRadius() sphere:
 \code
 const float zMin = zNearCoefficient() * zClippingCoefficient() * sceneRadius();
 if (zNear < zMin)
   zNear = zMin;
 // With an ORTHOGRAPHIC type, the value is simply clamped to 0.0
 \endcode

 See also the zFar(), zClippingCoefficient() and zNearCoefficient() documentations.

 If you need a completely different zNear computation, overload the zNear() and zFar() methods in a
 new class that publicly inherits from Camera and use QGLViewer::setCamera():
 \code
 class myCamera :: public qglviewer::Camera
 {
   virtual float Camera::zNear() const { return 0.001; };
   virtual float Camera::zFar() const { return 100.0; };
 }
 \endcode

 See the <a href="../examples/standardCamera.html">standardCamera example</a> for an application.

 \attention The value is always positive although the clipping plane is positioned at a negative z
 value in the Camera coordinate system. This follows the \c gluPerspective standard. */
float Camera::zNear() const
{
	const float zNearScene = zClippingCoefficient() * sceneRadius();
	float z = distanceToSceneCenter() - zNearScene;

	// Prevents negative or null zNear values.
	const float zMin = zNearCoefficient() * zNearScene;
	if (z < zMin)
		switch (type())
		{
		case Camera::PERSPECTIVE  : z = zMin; break;
		case Camera::ORTHOGRAPHIC : z = 0.0;  break;
		}
	return z;
}

/*! Returns the far clipping plane distance used by the Camera projection matrix.

The far clipping plane is positioned at a distance equal to zClippingCoefficient() * sceneRadius()
behind the sceneCenter():
\code
zFar = distanceToSceneCenter() + zClippingCoefficient()*sceneRadius();
\endcode

See the zNear() documentation for details. */
float Camera::zFar() const
{
	return distanceToSceneCenter() + zClippingCoefficient() * sceneRadius();
}


/*! Sets the vertical fieldOfView() of the Camera (in radians).

Note that focusDistance() is set to sceneRadius() / tan(fieldOfView()/2) by this method. */
void Camera::setFieldOfView(float fov) {	fieldOfView_ = fov;
	setFocusDistance(sceneRadius() / tan(fov/2.0));
	projectionMatrixIsUpToDate_ = false;
}

/*! Defines the Camera type().

Changing the camera Type alters the viewport and the objects' size can be changed. This method garantees that the two frustum match in a plane normal to viewDirection(), passing through the Revolve Around Point (RAP).

Prefix the type with \c Camera if needed, as in:
\code
camera()->setType(Camera::ORTHOGRAPHIC);
// or even qglviewer::Camera::ORTHOGRAPHIC if you do not use namespace
\endcode */
void Camera::setType(Type type)
{
	// make ORTHOGRAPHIC frustum fit PERSPECTIVE (at least in plane normal to viewDirection(), passing
	// through RAP). Done only when CHANGING type since orthoCoef_ may have been changed with a
	// setRevolveAroundPoint() in the meantime.
	if ( (type == Camera::ORTHOGRAPHIC) && (type_ == Camera::PERSPECTIVE) )
		orthoCoef_ = tan(fieldOfView()/2.0);
	type_ = type;
	projectionMatrixIsUpToDate_ = false;
}


void Camera::setFrame(ManipulatedCameraFrame* const mcf)
{
	if (!mcf)
		return;

	if (frame_) {
		disconnect(frame_, SIGNAL(modified()), this, SLOT(onFrameModified()));
	}

	frame_ = mcf;
	interpolationKfi_->setFrame(frame());

	connect(frame_, SIGNAL(modified()), this, SLOT(onFrameModified()));
	onFrameModified();
}

/*! Returns the distance from the Camera center to sceneCenter(), projected along the Camera Z axis.
  Used by zNear() and zFar() to optimize the Z range. */
float Camera::distanceToSceneCenter() const
{
	return fabs((frame()->coordinatesOf(sceneCenter())).z);
}


/*! Returns the \p halfWidth and \p halfHeight of the Camera orthographic frustum.

 These values are only valid and used when the Camera is of type() Camera::ORTHOGRAPHIC. They are
 expressed in OpenGL units and are used by loadProjectionMatrix() to define the projection matrix
 using:
 \code
 glOrtho( -halfWidth, halfWidth, -halfHeight, halfHeight, zNear(), zFar() )
 \endcode

 These values are proportional to the Camera (z projected) distance to the revolveAroundPoint().
 When zooming on the object, the Camera is translated forward \e and its frustum is narrowed, making
 the object appear bigger on screen, as intuitively expected.

 Overload this method to change this behavior if desired, as is done in the
 <a href="../examples/standardCamera.html">standardCamera example</a>. */
void Camera::getOrthoWidthHeight(GLdouble& halfWidth, GLdouble& halfHeight) const
{
	const float dist = orthoCoef_ * fabs(cameraCoordinatesOf(revolveAroundPoint()).z);
	//#CONNECTION# fitScreenRegion
	halfWidth  = dist * ((aspectRatio() < 1.0) ? 1.0 : aspectRatio());
	halfHeight = dist * ((aspectRatio() < 1.0) ? 1.0/aspectRatio() : 1.0);
}


/*! Computes the projection matrix associated with the Camera.

 If type() is Camera::PERSPECTIVE, defines a \c GL_PROJECTION matrix similar to what would \c
 gluPerspective() do using the fieldOfView(), window aspectRatio(), zNear() and zFar() parameters.

 If type() is Camera::ORTHOGRAPHIC, the projection matrix is as what \c glOrtho() would do.
 Frustum's width and height are set using getOrthoWidthHeight().

 Both types use zNear() and zFar() to place clipping planes. These values are determined from
 sceneRadius() and sceneCenter() so that they best fit the scene size.

 Use getProjectionMatrix() to retrieve this matrix. Overload loadProjectionMatrix() if you want your
 Camera to use an exotic projection matrix.

 \note You must call this method if your Camera is not associated with a QGLViewer and is used for
 offscreen computations (using (un)projectedCoordinatesOf() for instance). loadProjectionMatrix()
 does it otherwise. */
void Camera::computeProjectionMatrix() const
{
	if (projectionMatrixIsUpToDate_) return;

	const float ZNear = zNear();
	const float ZFar  = zFar();

	switch (type())
	{
	case Camera::PERSPECTIVE:
	{
		// #CONNECTION# all non null coefficients were set to 0.0 in constructor.
		const float f = 1.0/tan(fieldOfView()/2.0);
		projectionMatrix_[0]  = f/aspectRatio();
		projectionMatrix_[5]  = f;
		projectionMatrix_[10] = (ZNear + ZFar) / (ZNear - ZFar);
		projectionMatrix_[11] = -1.0;
		projectionMatrix_[14] = 2.0 * ZNear * ZFar / (ZNear - ZFar);
		projectionMatrix_[15] = 0.0;
		// same as gluPerspective( 180.0*fieldOfView()/M_PI, aspectRatio(), zNear(), zFar() );
		break;
	}
	case Camera::ORTHOGRAPHIC:
	{
		GLdouble w, h;
		getOrthoWidthHeight(w,h);
		projectionMatrix_[0]  = 1.0/w;
		projectionMatrix_[5]  = 1.0/h;
		projectionMatrix_[10] = -2.0/(ZFar - ZNear);
		projectionMatrix_[11] = 0.0;
		projectionMatrix_[14] = -(ZFar + ZNear)/(ZFar - ZNear);
		projectionMatrix_[15] = 1.0;
		// same as glOrtho( -w, w, -h, h, zNear(), zFar() );
		break;
	}
	}

	projectionMatrixIsUpToDate_ = true;
}

/*! Computes the modelView matrix associated with the Camera's position() and orientation().

 This matrix converts from the world coordinates system to the Camera coordinates system, so that
 coordinates can then be projected on screen using the projection matrix (see computeProjectionMatrix()).

 Use getModelViewMatrix() to retrieve this matrix.

 \note You must call this method if your Camera is not associated with a QGLViewer and is used for
 offscreen computations (using (un)projectedCoordinatesOf() for instance). loadModelViewMatrix()
 does it otherwise. */
void Camera::computeModelViewMatrix() const
{
	if (modelViewMatrixIsUpToDate_) return;

	const Quaternion q = frame()->orientation();

	const double q00 = 2.0l * q[0] * q[0];
	const double q11 = 2.0l * q[1] * q[1];
	const double q22 = 2.0l * q[2] * q[2];

	const double q01 = 2.0l * q[0] * q[1];
	const double q02 = 2.0l * q[0] * q[2];
	const double q03 = 2.0l * q[0] * q[3];

	const double q12 = 2.0l * q[1] * q[2];
	const double q13 = 2.0l * q[1] * q[3];

	const double q23 = 2.0l * q[2] * q[3];

	modelViewMatrix_[0] = 1.0l - q11 - q22;
	modelViewMatrix_[1] =        q01 - q23;
	modelViewMatrix_[2] =        q02 + q13;
	modelViewMatrix_[3] = 0.0l;

	modelViewMatrix_[4] =        q01 + q23;
	modelViewMatrix_[5] = 1.0l - q22 - q00;
	modelViewMatrix_[6] =        q12 - q03;
	modelViewMatrix_[7] = 0.0l;

	modelViewMatrix_[8] =        q02 - q13;
	modelViewMatrix_[9] =        q12 + q03;
	modelViewMatrix_[10] = 1.0l - q11 - q00;
	modelViewMatrix_[11] = 0.0l;

	const Vec t = q.inverseRotate(frame()->position());

	modelViewMatrix_[12] = -t.x;
	modelViewMatrix_[13] = -t.y;
	modelViewMatrix_[14] = -t.z;
	modelViewMatrix_[15] = 1.0l;

	modelViewMatrixIsUpToDate_ = true;
}


void Camera::loadProjectionMatrix(bool reset) const
{
	// WARNING: makeCurrent must be called by every calling method
	glMatrixMode(GL_PROJECTION);

	if (reset)
		glLoadIdentity();

	computeProjectionMatrix();

	glMultMatrixd(projectionMatrix_);
}


void Camera::loadModelViewMatrix(bool reset) const
{
	// WARNING: makeCurrent must be called by every calling method
	glMatrixMode(GL_MODELVIEW);
	computeModelViewMatrix();
	if (reset)
		glLoadMatrixd(modelViewMatrix_);
	else
		glMultMatrixd(modelViewMatrix_);
}

void Camera::loadProjectionMatrixStereo(bool leftBuffer) const
{
	float left, right, bottom, top;
	float screenHalfWidth, halfWidth, side, shift, delta;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	switch (type())
	{
	case Camera::PERSPECTIVE:
		// compute half width of screen,
		// corresponding to zero parallax plane to deduce decay of cameras
		screenHalfWidth = focusDistance() * tan(horizontalFieldOfView() / 2.0);
		shift = screenHalfWidth * IODistance() / physicalScreenWidth();
		// should be * current y  / y total
		// to take into account that the window doesn't cover the entire screen

		// compute half width of "view" at znear and the delta corresponding to
		// the shifted camera to deduce what to set for asymmetric frustums
		halfWidth = zNear() * tan(horizontalFieldOfView() / 2.0);
		delta  = shift * zNear() / focusDistance();
		side   = leftBuffer ? -1.0 : 1.0;

		left   = -halfWidth + side * delta;
		right  =  halfWidth + side * delta;
		top    = halfWidth / aspectRatio();
		bottom = -top;
		glFrustum(left, right, bottom, top, zNear(), zFar() );
		break;

	case Camera::ORTHOGRAPHIC:
		qWarning("Camera::setProjectionMatrixStereo: Stereo not available with Ortho mode");
		break;
	}
}

void Camera::loadModelViewMatrixStereo(bool leftBuffer) const
{
	// WARNING: makeCurrent must be called by every calling method
	glMatrixMode(GL_MODELVIEW);

	float halfWidth = focusDistance() * tan(horizontalFieldOfView() / 2.0);
	float shift     = halfWidth * IODistance() / physicalScreenWidth(); // * current window width / full screen width

	computeModelViewMatrix();
	if (leftBuffer)
		modelViewMatrix_[12] -= shift;
	else
		modelViewMatrix_[12] += shift;
	glLoadMatrixd(modelViewMatrix_);
}

void Camera::getProjectionMatrix(GLdouble m[16]) const
{
	computeProjectionMatrix();
	for (unsigned short i=0; i<16; ++i)
		m[i] = projectionMatrix_[i];
}


void Camera::getModelViewMatrix(GLdouble m[16]) const
{
	// May not be needed, but easier like this.
	// Prevents from retrieving matrix in stereo mode -> overwrites shifted value.
	computeModelViewMatrix();
	for (unsigned short i=0; i<16; ++i)
		m[i] = modelViewMatrix_[i];
}

void Camera::getModelViewProjectionMatrix(GLdouble m[16]) const
{
	GLdouble mv[16];
	GLdouble proj[16];
	getModelViewMatrix(mv);
	getProjectionMatrix(proj);

	for (unsigned short i=0; i<4; ++i)
	{
		for (unsigned short j=0; j<4; ++j)
		{
			double sum = 0.0;
			for (unsigned short k=0; k<4; ++k)
				sum += proj[i+4*k]*mv[k+4*j];
			m[i+4*j] = sum;
		}
	}
}

#ifndef DOXYGEN
void Camera::getProjectionMatrix(GLfloat m[16]) const
{
	qWarning("Warning : Camera::getProjectionMatrix requires a GLdouble matrix array");
	static GLdouble mat[16];
	getProjectionMatrix(mat);
	for (int i=0; i<16; ++i)
		m[i] = float(mat[i]);
}

void Camera::getModelViewMatrix(GLfloat m[16]) const
{
	qWarning("Warning : Camera::getModelViewMatrix requires a GLdouble matrix array");
	static GLdouble mat[16];
	getModelViewMatrix(mat);
	for (int i=0; i<16; ++i)
		m[i] = float(mat[i]);
}
#endif


void Camera::setSceneRadius(float radius)
{
	if (radius <= 0.0)
	{
		qWarning("Scene radius must be positive - Ignoring value");
		return;
	}

	sceneRadius_ = radius;
	projectionMatrixIsUpToDate_ = false;

	setFocusDistance(sceneRadius() / tan(fieldOfView()/2.0));

	frame()->setFlySpeed(0.01*sceneRadius());
}

/*! Similar to setSceneRadius() and setSceneCenter(), but the scene limits are defined by a (world
  axis aligned) bounding box. */
void Camera::setSceneBoundingBox(const Vec& min, const Vec& max)
{
	setSceneCenter((min+max)/2.0);
	setSceneRadius(0.5*(max-min).norm());
}


/*! Sets the sceneCenter().

 \attention This method also sets the revolveAroundPoint() to sceneCenter(). */
void Camera::setSceneCenter(const Vec& center)
{
	sceneCenter_ = center;
	setRevolveAroundPoint(sceneCenter());
	projectionMatrixIsUpToDate_ = false;
}

bool Camera::setSceneCenterFromPixel(const QPoint& pixel)
{
	bool found;
	Vec point = pointUnderPixel(pixel, found);
	if (found)
		setSceneCenter(point);
	return found;
}

void Camera::setRevolveAroundPoint(const Vec& rap)
{
	const float prevDist = fabs(cameraCoordinatesOf(revolveAroundPoint()).z);

	// If frame's RAP is set directly, projectionMatrixIsUpToDate_ should also be
	// set to false to ensure proper recomputation of the ORTHO projection matrix.
	frame()->setRevolveAroundPoint(rap);

	// orthoCoef_ is used to compensate for changes of the revolveAroundPoint, so that the image does
	// not change when the revolveAroundPoint is changed in ORTHOGRAPHIC mode.
	const float newDist = fabs(cameraCoordinatesOf(revolveAroundPoint()).z);
	// Prevents division by zero when rap is set to camera position
	if ((prevDist > 1E-9) && (newDist > 1E-9))
		orthoCoef_ *= prevDist / newDist;
	projectionMatrixIsUpToDate_ = false;
}

bool Camera::setRevolveAroundPointFromPixel(const QPoint& pixel)
{
	bool found;
	Vec point = pointUnderPixel(pixel, found);
	if (found)
		setRevolveAroundPoint(point);
	return found;
}

float Camera::pixelGLRatio(const Vec& position) const
{
	switch (type())
	{
	case Camera::PERSPECTIVE :
		return 2.0 * fabs((frame()->coordinatesOf(position)).z) * tan(fieldOfView()/2.0) / screenHeight();
	case Camera::ORTHOGRAPHIC :
	{
		GLdouble w, h;
		getOrthoWidthHeight(w,h);
		return 2.0 * h / screenHeight();
	}
	}
	// Bad compilers complain
	return 1.0;
}

void Camera::setFOVToFitScene()
{
	if (distanceToSceneCenter() > sqrt(2.0)*sceneRadius())
		setFieldOfView(2.0 * asin(sceneRadius() / distanceToSceneCenter()));
	else
		setFieldOfView(M_PI / 2.0f);
}


void Camera::interpolateToZoomOnPixel(const QPoint& pixel)
{
	const float coef = 0.1f;

	bool found;
	Vec target = pointUnderPixel(pixel, found);

	if (!found)
		return;

	if (interpolationKfi_->interpolationIsStarted())
		interpolationKfi_->stopInterpolation();

	interpolationKfi_->deletePath();
	interpolationKfi_->addKeyFrame(*(frame()));

	interpolationKfi_->addKeyFrame(Frame(0.3f*frame()->position() + 0.7f*target, frame()->orientation()), 0.4f);

	// Small hack: attach a temporary frame to take advantage of lookAt without modifying frame
	static ManipulatedCameraFrame* tempFrame = new ManipulatedCameraFrame();
	ManipulatedCameraFrame* const originalFrame = frame();
	tempFrame->setPosition(coef*frame()->position() + (1.0-coef)*target);
	tempFrame->setOrientation(frame()->orientation());
	setFrame(tempFrame);
	lookAt(target);
	setFrame(originalFrame);

	interpolationKfi_->addKeyFrame(*(tempFrame), 1.0);

	interpolationKfi_->startInterpolation();
}


void Camera::interpolateToFitScene()
{
	if (interpolationKfi_->interpolationIsStarted())
		interpolationKfi_->stopInterpolation();

	interpolationKfi_->deletePath();
	interpolationKfi_->addKeyFrame(*(frame()));

	// Small hack:  attach a temporary frame to take advantage of lookAt without modifying frame
	static ManipulatedCameraFrame* tempFrame = new ManipulatedCameraFrame();
	ManipulatedCameraFrame* const originalFrame = frame();
	tempFrame->setPosition(frame()->position());
	tempFrame->setOrientation(frame()->orientation());
	setFrame(tempFrame);
	showEntireScene();
	setFrame(originalFrame);

	interpolationKfi_->addKeyFrame(*(tempFrame));

	interpolationKfi_->startInterpolation();
}


void Camera::interpolateTo(const Frame& fr, float duration)
{
	if (interpolationKfi_->interpolationIsStarted())
		interpolationKfi_->stopInterpolation();

	interpolationKfi_->deletePath();
	interpolationKfi_->addKeyFrame(*(frame()));
	interpolationKfi_->addKeyFrame(fr, duration);

	interpolationKfi_->startInterpolation();
}


Vec Camera::pointUnderPixel(const QPoint& pixel, bool& found) const
{
	float depth;
	// Qt uses upper corner for its origin while GL uses the lower corner.
	glReadPixels(pixel.x(), screenHeight()-1-pixel.y(), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
	found = depth < 1.0;
	Vec point(pixel.x(), pixel.y(), depth);
	point = unprojectedCoordinatesOf(point);
	return point;
}

void Camera::showEntireScene()
{
	fitSphere(sceneCenter(), sceneRadius());
}


void Camera::centerScene()
{
	frame()->projectOnLine(sceneCenter(), viewDirection());
}


void Camera::lookAt(const Vec& target)
{
	setViewDirection(target - position());
}


void Camera::fitSphere(const Vec& center, float radius)
{
	float distance = 0.0f;
	switch (type())
	{
	case Camera::PERSPECTIVE :
	{
		const float yview = radius / sin(fieldOfView()/2.0);
		const float xview = radius / sin(horizontalFieldOfView()/2.0);
		distance = qMax(xview,yview);
		break;
	}
	case Camera::ORTHOGRAPHIC :
	{
		distance = ((center-revolveAroundPoint()) * viewDirection()) + (radius / orthoCoef_);
		break;
	}
	}
	Vec newPos(center - distance * viewDirection());
	frame()->setPositionWithConstraint(newPos);
}


void Camera::fitBoundingBox(const Vec& min, const Vec& max)
{
	float diameter = qMax(fabs(max[1]-min[1]), fabs(max[0]-min[0]));
	diameter = qMax(fabsf(max[2]-min[2]), diameter);
	fitSphere(0.5*(min+max), 0.5*diameter);
}


void Camera::fitScreenRegion(const QRect& rectangle)
{
	const Vec vd = viewDirection();
	const float distToPlane = distanceToSceneCenter();
	const QPoint center = rectangle.center();

	Vec orig, dir;
	convertClickToLine( center, orig, dir );
	Vec newCenter = orig + distToPlane / (dir*vd) * dir;

	convertClickToLine( QPoint(rectangle.x(), center.y()), orig, dir );
	const Vec pointX = orig + distToPlane / (dir*vd) * dir;

	convertClickToLine( QPoint(center.x(), rectangle.y()), orig, dir );
	const Vec pointY = orig + distToPlane / (dir*vd) * dir;

	float distance = 0.0f;
	switch (type())
	{
	case Camera::PERSPECTIVE :
	{
		const float distX = (pointX-newCenter).norm() / sin(horizontalFieldOfView()/2.0);
		const float distY = (pointY-newCenter).norm() / sin(fieldOfView()/2.0);
		distance = qMax(distX, distY);
		break;
	}
	case Camera::ORTHOGRAPHIC :
	{
		const float dist = ((newCenter-revolveAroundPoint()) * vd);
		//#CONNECTION# getOrthoWidthHeight
		const float distX = (pointX-newCenter).norm() / orthoCoef_ / ((aspectRatio() < 1.0) ? 1.0 : aspectRatio());
		const float distY = (pointY-newCenter).norm() / orthoCoef_ / ((aspectRatio() < 1.0) ? 1.0/aspectRatio() : 1.0);
		distance = dist + qMax(distX, distY);
		break;
	}
	}

	Vec newPos(newCenter - distance * vd);
	frame()->setPositionWithConstraint(newPos);
}


void Camera::setUpVector(const Vec& up, bool noMove)
{
	Quaternion q(Vec(0.0, 1.0, 0.0), frame()->transformOf(up));

	if (!noMove)
		frame()->setPosition(revolveAroundPoint() - (frame()->orientation()*q).rotate(frame()->coordinatesOf(revolveAroundPoint())));

	frame()->rotate(q);

	// Useful in fly mode to keep the horizontal direction.
	frame()->updateFlyUpVector();
}

void Camera::setOrientation(float theta, float phi)
{
	Vec axis(0.0, 1.0, 0.0);
	const Quaternion rot1(axis, theta);
	axis = Vec(-cos(theta), 0., sin(theta));
	const Quaternion rot2(axis, phi);
	setOrientation(rot1 * rot2);
}

/*! Sets the Camera orientation(), defined in the world coordinate system. */
void Camera::setOrientation(const Quaternion& q)
{
	frame()->setOrientation(q);
	frame()->updateFlyUpVector();
}

void Camera::setViewDirection(const Vec& direction)
{
	if (direction.squaredNorm() < 1E-10)
		return;

	Vec xAxis = direction ^ upVector();
	if (xAxis.squaredNorm() < 1E-10)
	{
		// target is aligned with upVector, this means a rotation around X axis
		// X axis is then unchanged, let's keep it !
		xAxis = frame()->inverseTransformOf(Vec(1.0, 0.0, 0.0));
	}

	Quaternion q;
	q.setFromRotatedBasis(xAxis, xAxis^direction, -direction);
	frame()->setOrientationWithConstraint(q);
}

// Compute a 3 by 3 determinant.
static float det(float m00,float m01,float m02,
				 float m10,float m11,float m12,
				 float m20,float m21,float m22)
{
	return m00*m11*m22 + m01*m12*m20 + m02*m10*m21 - m20*m11*m02 - m10*m01*m22 - m00*m21*m12;
}

// Computes the index of element [i][j] in a \c float matrix[3][4].
static inline unsigned int ind(unsigned int i, unsigned int j)
{
	return (i*4+j);
}


Vec Camera::position() const { return frame()->position(); }


Vec Camera::upVector() const
{
	return frame()->inverseTransformOf(Vec(0.0, 1.0, 0.0));
}

Vec Camera::viewDirection() const { return frame()->inverseTransformOf(Vec(0.0, 0.0, -1.0)); }


Vec Camera::rightVector() const
{
	return frame()->inverseTransformOf(Vec(1.0, 0.0, 0.0));
}

Quaternion Camera::orientation() const { return frame()->orientation(); }


void Camera::setPosition(const Vec& pos) { frame()->setPosition(pos); }

/*! Returns the Camera frame coordinates of a point \p src defined in world coordinates.

worldCoordinatesOf() performs the inverse transformation.

Note that the point coordinates are simply converted in a different coordinate system. They are
not projected on screen. Use projectedCoordinatesOf() for that. */
Vec Camera::cameraCoordinatesOf(const Vec& src) const { return frame()->coordinatesOf(src); }

/*! Returns the world coordinates of the point whose position \p src is defined in the Camera
coordinate system.

cameraCoordinatesOf() performs the inverse transformation. */
Vec Camera::worldCoordinatesOf(const Vec& src) const { return frame()->inverseCoordinatesOf(src); }

/*! Returns the fly speed of the Camera.

Simply returns frame()->flySpeed(). See the ManipulatedCameraFrame::flySpeed() documentation.
This value is only meaningful when the MouseAction bindings is QGLViewer::MOVE_FORWARD or
QGLViewer::MOVE_BACKWARD.

Set to 1% of the sceneRadius() by setSceneRadius(). See also setFlySpeed(). */
float Camera::flySpeed() const { return frame()->flySpeed(); }

/*! Sets the Camera flySpeed().

\attention This value is modified by setSceneRadius(). */
void Camera::setFlySpeed(float speed) { frame()->setFlySpeed(speed); }

/*! The point the Camera revolves around with the QGLViewer::ROTATE mouse binding. Defined in world coordinate system.

Default value is the sceneCenter().

\attention setSceneCenter() changes this value. */
Vec Camera::revolveAroundPoint() const { return frame()->revolveAroundPoint(); }


void Camera::setFromModelViewMatrix(const GLdouble* const modelViewMatrix)
{
	// Get upper left (rotation) matrix
	double upperLeft[3][3];
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			upperLeft[i][j] = modelViewMatrix[i*4+j];

	// Transform upperLeft into the associated Quaternion
	Quaternion q;
	q.setFromRotationMatrix(upperLeft);

	setOrientation(q);
	setPosition(-q.rotate(Vec(modelViewMatrix[12], modelViewMatrix[13], modelViewMatrix[14])));
}


void Camera::setFromProjectionMatrix(const float matrix[12])
{
	// The 3 lines of the matrix are the normals to the planes x=0, y=0, z=0
	// in the camera CS. As we normalize them, we do not need the 4th coordinate.
	Vec line_0(matrix[ind(0,0)],matrix[ind(0,1)],matrix[ind(0,2)]);
	Vec line_1(matrix[ind(1,0)],matrix[ind(1,1)],matrix[ind(1,2)]);
	Vec line_2(matrix[ind(2,0)],matrix[ind(2,1)],matrix[ind(2,2)]);

	line_0.normalize();
	line_1.normalize();
	line_2.normalize();

	const Vec cam_pos = Vec(det(matrix[ind(0,1)],matrix[ind(0,2)],matrix[ind(0,3)],
			matrix[ind(1,1)],matrix[ind(1,2)],matrix[ind(1,3)],
			matrix[ind(2,1)],matrix[ind(2,2)],matrix[ind(2,3)]),

			-det(matrix[ind(0,0)],matrix[ind(0,2)],matrix[ind(0,3)],
			matrix[ind(1,0)],matrix[ind(1,2)],matrix[ind(1,3)],
			matrix[ind(2,0)],matrix[ind(2,2)],matrix[ind(2,3)]),

			det(matrix[ind(0,0)],matrix[ind(0,1)],matrix[ind(0,3)],
			matrix[ind(1,0)],matrix[ind(1,1)],matrix[ind(1,3)],
			matrix[ind(2,0)],matrix[ind(2,1)],matrix[ind(2,3)])) /

			(-det(matrix[ind(0,0)],matrix[ind(0,1)],matrix[ind(0,2)],
			matrix[ind(1,0)],matrix[ind(1,1)],matrix[ind(1,2)],
			matrix[ind(2,0)],matrix[ind(2,1)],matrix[ind(2,2)]));

	// We compute the rotation matrix column by column.

	// GL Z axis is front facing.
	Vec column_2 = -line_2;

	// X-axis is almost like line_0 but should be orthogonal to the Z axis.
	Vec column_0 = ((column_2^line_0)^column_2);
	column_0.normalize();

	// Y-axis is almost like line_1 but should be orthogonal to the Z axis.
	// Moreover line_1 is downward oriented as the screen CS.
	Vec column_1 = -((column_2^line_1)^column_2);
	column_1.normalize();

	double rot[3][3];
	rot[0][0] = column_0[0];
	rot[1][0] = column_0[1];
	rot[2][0] = column_0[2];

	rot[0][1] = column_1[0];
	rot[1][1] = column_1[1];
	rot[2][1] = column_1[2];

	rot[0][2] = column_2[0];
	rot[1][2] = column_2[1];
	rot[2][2] = column_2[2];

	// We compute the field of view

	// line_1^column_0 -> vector of intersection line between
	// y_screen=0 and x_camera=0 plane.
	// column_2*(...)  -> cos of the angle between Z vector et y_screen=0 plane
	// * 2 -> field of view = 2 * half angle

	// We need some intermediate values.
	Vec dummy = line_1^column_0;
	dummy.normalize();
	float fov = acos(column_2*dummy) * 2.0;

	// We set the camera.
	Quaternion q;
	q.setFromRotationMatrix(rot);
	setOrientation(q);
	setPosition(cam_pos);
	setFieldOfView(fov);
}



void Camera::getCameraCoordinatesOf(const float src[3], float res[3]) const
{
	Vec r = cameraCoordinatesOf(Vec(src));
	for (int i=0; i<3; ++i)
		res[i] = r[i];
}

/*! Same as worldCoordinatesOf(), but with \c float[3] parameters (\p src and \p res may be identical pointers). */
void Camera::getWorldCoordinatesOf(const float src[3], float res[3]) const
{
	Vec r = worldCoordinatesOf(Vec(src));
	for (int i=0; i<3; ++i)
		res[i] = r[i];
}

/*! Fills \p viewport with the Camera OpenGL viewport.

This method is mainly used in conjunction with \c gluProject, which requires such a viewport.
Returned values are (0, screenHeight(), screenWidth(), - screenHeight()), so that the origin is
located in the \e upper left corner of the window (Qt style coordinate system). */
void Camera::getViewport(GLint viewport[4]) const
{
	viewport[0] = 0;
	viewport[1] = screenHeight();
	viewport[2] = screenWidth();
	viewport[3] = -screenHeight();
}


Vec Camera::projectedCoordinatesOf(const Vec& src, const Frame* frame) const
{
	GLdouble x,y,z;
	static GLint viewport[4];
	getViewport(viewport);

	if (frame)
	{
		const Vec tmp = frame->inverseCoordinatesOf(src);
		gluProject(tmp.x,tmp.y,tmp.z, modelViewMatrix_, projectionMatrix_, viewport,  &x,&y,&z);
	}
	else
		gluProject(src.x,src.y,src.z, modelViewMatrix_, projectionMatrix_, viewport,  &x,&y,&z);

	return Vec(x,y,z);
}
Vec Camera::unprojectedCoordinatesOf(const Vec& src, const Frame* frame) const
{
	GLdouble x,y,z;
	static GLint viewport[4];
	getViewport(viewport);
	gluUnProject(src.x,src.y,src.z, modelViewMatrix_,  projectionMatrix_,  viewport,  &x,&y,&z);
	if (frame)
		return frame->coordinatesOf(Vec(x,y,z));
	else
		return Vec(x,y,z);
}

/*! Same as projectedCoordinatesOf(), but with \c float parameters (\p src and \p res can be identical pointers). */
void Camera::getProjectedCoordinatesOf(const float src[3], float res[3], const Frame* frame) const
{
	Vec r = projectedCoordinatesOf(Vec(src), frame);
	for (int i=0; i<3; ++i)
		res[i] = r[i];
}

/*! Same as unprojectedCoordinatesOf(), but with \c float parameters (\p src and \p res can be identical pointers). */
void Camera::getUnprojectedCoordinatesOf(const float src[3], float res[3], const Frame* frame) const
{
	Vec r = unprojectedCoordinatesOf(Vec(src), frame);
	for (int i=0; i<3; ++i)
		res[i] = r[i];
}


KeyFrameInterpolator* Camera::keyFrameInterpolator(int i) const
{
	if (kfi_.contains(i))
		return kfi_[i];
	else
		return NULL;
}


void Camera::setKeyFrameInterpolator(int i, KeyFrameInterpolator* const kfi)
{
	if (kfi)
		kfi_[i] = kfi;
	else
		kfi_.remove(i);
}


void Camera::addKeyFrameToPath(int i)
{
	if (!kfi_.contains(i))
		setKeyFrameInterpolator(i, new KeyFrameInterpolator(frame()));

	kfi_[i]->addKeyFrame(*(frame()));
}


void Camera::playPath(int i)
{
	if (kfi_.contains(i)) {
		if (kfi_[i]->interpolationIsStarted())
			kfi_[i]->stopInterpolation();
		else
			kfi_[i]->startInterpolation();
	}
}


void Camera::resetPath(int i)
{
	if (kfi_.contains(i)) {
		if ((kfi_[i]->interpolationIsStarted()))
			kfi_[i]->stopInterpolation();
		else
		{
			kfi_[i]->resetInterpolation();
			kfi_[i]->interpolateAtTime(kfi_[i]->interpolationTime());
		}
	}
}


void Camera::deletePath(int i)
{
	if (kfi_.contains(i))
	{
		kfi_[i]->stopInterpolation();
		delete kfi_[i];
		kfi_.remove(i);
	}
}


void Camera::drawAllPaths()
{
	for (QMap<int, KeyFrameInterpolator*>::ConstIterator it = kfi_.begin(), end=kfi_.end(); it != end; ++it)
		(it.value())->drawPath(3, 5, sceneRadius());
}


QDomElement Camera::domElement(const QString& name, QDomDocument& document) const
{
	QDomElement de = document.createElement(name);
	QDomElement paramNode = document.createElement("Parameters");
	paramNode.setAttribute("fieldOfView", QString::number(fieldOfView()));
	paramNode.setAttribute("zNearCoefficient", QString::number(zNearCoefficient()));
	paramNode.setAttribute("zClippingCoefficient", QString::number(zClippingCoefficient()));
	paramNode.setAttribute("orthoCoef", QString::number(orthoCoef_));
	paramNode.setAttribute("sceneRadius", QString::number(sceneRadius()));
	paramNode.appendChild(sceneCenter().domElement("SceneCenter", document));

	switch (type())
	{
	case Camera::PERSPECTIVE  :	paramNode.setAttribute("Type", "PERSPECTIVE"); break;
	case Camera::ORTHOGRAPHIC :	paramNode.setAttribute("Type", "ORTHOGRAPHIC"); break;
	}
	de.appendChild(paramNode);

	QDomElement stereoNode = document.createElement("Stereo");
	stereoNode.setAttribute("IODist", QString::number(IODistance()));
	stereoNode.setAttribute("focusDistance", QString::number(focusDistance()));
	stereoNode.setAttribute("physScreenWidth", QString::number(physicalScreenWidth()));
	de.appendChild(stereoNode);

	de.appendChild(frame()->domElement("ManipulatedCameraFrame", document));

	// KeyFrame paths
	for (QMap<int, KeyFrameInterpolator*>::ConstIterator it = kfi_.begin(), end=kfi_.end(); it != end; ++it)
	{
		QDomElement kfNode = (it.value())->domElement("KeyFrameInterpolator", document);
		kfNode.setAttribute("index", QString::number(it.key()));
		de.appendChild(kfNode);
	}

	return de;
}


void Camera::initFromDOMElement(const QDomElement& element)
{
	QDomElement child=element.firstChild().toElement();

	QMutableMapIterator<int, KeyFrameInterpolator*> it(kfi_);
	while (it.hasNext()) {
		it.next();
		deletePath(it.key());
	}

	while (!child.isNull())
	{
		if (child.tagName() == "Parameters")
		{
			// #CONNECTION# Default values set in constructor
			setFieldOfView(DomUtils::floatFromDom(child, "fieldOfView", M_PI/4.0f));
			setZNearCoefficient(DomUtils::floatFromDom(child, "zNearCoefficient", 0.005f));
			setZClippingCoefficient(DomUtils::floatFromDom(child, "zClippingCoefficient", sqrt(3.0)));
			orthoCoef_ = DomUtils::floatFromDom(child, "orthoCoef", tan(fieldOfView()/2.0));
			setSceneRadius(DomUtils::floatFromDom(child, "sceneRadius", sceneRadius()));

			setType(PERSPECTIVE);
			QString type = child.attribute("Type", "PERSPECTIVE");
			if (type == "PERSPECTIVE")  setType(Camera::PERSPECTIVE);
			if (type == "ORTHOGRAPHIC") setType(Camera::ORTHOGRAPHIC);

			QDomElement child2=child.firstChild().toElement();
			while (!child2.isNull())
			{
				/* Although the scene does not change when a camera is loaded, restore the saved center and radius values.
		   Mainly useful when a the viewer is restored on startup, with possible additional cameras. */
				if (child2.tagName() == "SceneCenter")
					setSceneCenter(Vec(child2));

				child2 = child2.nextSibling().toElement();
			}
		}

		if (child.tagName() == "ManipulatedCameraFrame")
			frame()->initFromDOMElement(child);

		if (child.tagName() == "Stereo")
		{
			setIODistance(DomUtils::floatFromDom(child, "IODist", 0.062f));
			setFocusDistance(DomUtils::floatFromDom(child, "focusDistance", focusDistance()));
			setPhysicalScreenWidth(DomUtils::floatFromDom(child, "physScreenWidth", 0.5f));
		}

		if (child.tagName() == "KeyFrameInterpolator")
		{
			int index = DomUtils::intFromDom(child, "index", 0);
			setKeyFrameInterpolator(index, new KeyFrameInterpolator(frame()));
			if (keyFrameInterpolator(index))
				keyFrameInterpolator(index)->initFromDOMElement(child);
		}

		child = child.nextSibling().toElement();
	}
}


void Camera::convertClickToLine(const QPoint& pixel, Vec& orig, Vec& dir) const
{
	switch (type())
	{
	case Camera::PERSPECTIVE:
		orig = position();
		dir = Vec( ((2.0 * pixel.x() / screenWidth()) - 1.0) * tan(fieldOfView()/2.0) * aspectRatio(),
				   ((2.0 * (screenHeight()-pixel.y()) / screenHeight()) - 1.0) * tan(fieldOfView()/2.0),
				   -1.0 );
		dir = worldCoordinatesOf(dir) - orig;
		dir.normalize();
		break;

	case Camera::ORTHOGRAPHIC:
	{
		GLdouble w,h;
		getOrthoWidthHeight(w,h);
		orig = Vec((2.0 * pixel.x() / screenWidth() - 1.0)*w, -(2.0 * pixel.y() / screenHeight() - 1.0)*h, 0.0);
		orig = worldCoordinatesOf(orig);
		dir = viewDirection();
		break;
	}
	}
}

#ifndef DOXYGEN
/*! This method has been deprecated in libQGLViewer version 2.2.0 */
void Camera::drawCamera(float, float, float)
{
	qWarning("drawCamera is deprecated. Use Camera::draw() instead.");
}
#endif


void Camera::draw(bool drawFarPlane, float scale) const
{
	glPushMatrix();
	glMultMatrixd(frame()->worldMatrix());

	// 0 is the upper left coordinates of the near corner, 1 for the far one
	Vec points[2];

	points[0].z = scale * zNear();
	points[1].z = scale * zFar();

	switch (type())
	{
	case Camera::PERSPECTIVE:
	{
		points[0].y = points[0].z * tan(fieldOfView()/2.0);
		points[0].x = points[0].y * aspectRatio();

		const float ratio = points[1].z / points[0].z;

		points[1].y = ratio * points[0].y;
		points[1].x = ratio * points[0].x;
		break;
	}
	case Camera::ORTHOGRAPHIC:
	{
		GLdouble hw, hh;
		getOrthoWidthHeight(hw, hh);
		points[0].x = points[1].x = scale * float(hw);
		points[0].y = points[1].y = scale * float(hh);
		break;
	}
	}

	const int farIndex = drawFarPlane?1:0;

	// Near and (optionally) far plane(s)
	glBegin(GL_QUADS);
	for (int i=farIndex; i>=0; --i)
	{
		glNormal3f(0.0, 0.0, (i==0)?1.0:-1.0);
		glVertex3f( points[i].x,  points[i].y, -points[i].z);
		glVertex3f(-points[i].x,  points[i].y, -points[i].z);
		glVertex3f(-points[i].x, -points[i].y, -points[i].z);
		glVertex3f( points[i].x, -points[i].y, -points[i].z);
	}
	glEnd();

	// Up arrow
	const float arrowHeight    = 1.5f * points[0].y;
	const float baseHeight     = 1.2f * points[0].y;
	const float arrowHalfWidth = 0.5f * points[0].x;
	const float baseHalfWidth  = 0.3f * points[0].x;

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	// Base
	glBegin(GL_QUADS);
	glVertex3f(-baseHalfWidth, points[0].y, -points[0].z);
	glVertex3f( baseHalfWidth, points[0].y, -points[0].z);
	glVertex3f( baseHalfWidth, baseHeight,  -points[0].z);
	glVertex3f(-baseHalfWidth, baseHeight,  -points[0].z);
	glEnd();

	// Arrow
	glBegin(GL_TRIANGLES);
	glVertex3f( 0.0f,           arrowHeight, -points[0].z);
	glVertex3f(-arrowHalfWidth, baseHeight,  -points[0].z);
	glVertex3f( arrowHalfWidth, baseHeight,  -points[0].z);
	glEnd();

	// Frustum lines
	switch (type())
	{
	case Camera::PERSPECTIVE :
		glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f( points[farIndex].x,  points[farIndex].y, -points[farIndex].z);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(-points[farIndex].x,  points[farIndex].y, -points[farIndex].z);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(-points[farIndex].x, -points[farIndex].y, -points[farIndex].z);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f( points[farIndex].x, -points[farIndex].y, -points[farIndex].z);
		glEnd();
		break;
	case Camera::ORTHOGRAPHIC :
		if (drawFarPlane)
		{
			glBegin(GL_LINES);
			glVertex3f( points[0].x,  points[0].y, -points[0].z);
			glVertex3f( points[1].x,  points[1].y, -points[1].z);
			glVertex3f(-points[0].x,  points[0].y, -points[0].z);
			glVertex3f(-points[1].x,  points[1].y, -points[1].z);
			glVertex3f(-points[0].x, -points[0].y, -points[0].z);
			glVertex3f(-points[1].x, -points[1].y, -points[1].z);
			glVertex3f( points[0].x, -points[0].y, -points[0].z);
			glVertex3f( points[1].x, -points[1].y, -points[1].z);
			glEnd();
		}
	}

	glPopMatrix();
}


void Camera::getFrustumPlanesCoefficients(GLdouble coef[6][4]) const
{
	// Computed once and for all
	const Vec pos          = position();
	const Vec viewDir      = viewDirection();
	const Vec up           = upVector();
	const Vec right        = rightVector();
	const float posViewDir = pos * viewDir;

	static Vec normal[6];
	static GLdouble dist[6];

	switch (type())
	{
	case Camera::PERSPECTIVE :
	{
		const float hhfov = horizontalFieldOfView() / 2.0;
		const float chhfov = cos(hhfov);
		const float shhfov = sin(hhfov);
		normal[0] = - shhfov * viewDir;
		normal[1] = normal[0] + chhfov * right;
		normal[0] = normal[0] - chhfov * right;

		normal[2] = -viewDir;
		normal[3] =  viewDir;

		const float hfov = fieldOfView() / 2.0;
		const float chfov = cos(hfov);
		const float shfov = sin(hfov);
		normal[4] = - shfov * viewDir;
		normal[5] = normal[4] - chfov * up;
		normal[4] = normal[4] + chfov * up;

		for (int i=0; i<2; ++i)
			dist[i] = pos * normal[i];
		for (int j=4; j<6; ++j)
			dist[j] = pos * normal[j];

		// Natural equations are:
		// dist[0,1,4,5] = pos * normal[0,1,4,5];
		// dist[2] = (pos + zNear() * viewDir) * normal[2];
		// dist[3] = (pos + zFar()  * viewDir) * normal[3];

		// 2 times less computations using expanded/merged equations. Dir vectors are normalized.
		const float posRightCosHH = chhfov * pos * right;
		dist[0] = -shhfov * posViewDir;
		dist[1] = dist[0] + posRightCosHH;
		dist[0] = dist[0] - posRightCosHH;
		const float posUpCosH = chfov * pos * up;
		dist[4] = - shfov * posViewDir;
		dist[5] = dist[4] - posUpCosH;
		dist[4] = dist[4] + posUpCosH;

		break;
	}
	case Camera::ORTHOGRAPHIC :
		normal[0] = -right;
		normal[1] =  right;
		normal[4] =  up;
		normal[5] = -up;

		GLdouble hw, hh;
		getOrthoWidthHeight(hw, hh);
		dist[0] = (pos - hw * right) * normal[0];
		dist[1] = (pos + hw * right) * normal[1];
		dist[4] = (pos + hh * up) * normal[4];
		dist[5] = (pos - hh * up) * normal[5];
		break;
	}

	// Front and far planes are identical for both camera types.
	normal[2] = -viewDir;
	normal[3] =  viewDir;
	dist[2] = -posViewDir - zNear();
	dist[3] =  posViewDir + zFar();

	for (int i=0; i<6; ++i)
	{
		coef[i][0] = GLdouble(normal[i].x);
		coef[i][1] = GLdouble(normal[i].y);
		coef[i][2] = GLdouble(normal[i].z);
		coef[i][3] = dist[i];
	}
}

void Camera::onFrameModified() {
	projectionMatrixIsUpToDate_ = false;
	modelViewMatrixIsUpToDate_ = false;
}
