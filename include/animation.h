#ifndef ANIMATION_H
#define ANIMATION_H

#include <skeleton.h>
#include <rbf.h>

/**
 * \class Animation
 *
 * \brief A handler for the animation process
 *
 * This class encapsulates the basic animation process. The animation 
 * process in short is: record a skeleton frame, associate it with a 
 * marker position and later move the controller between the markers to
 * get an interpolated pose. The animation builds the interpolated 
 * skeleton.
 *
 * \note current implementation is with a single controller.
 */

class Animation {
    protected:
    /// current interpolated skeleton
    Skeleton skel;
    /// array of markers
    std::vector<Eigen::Vector4d> markers;
    /// controller for interpolating between markers
    Eigen::Vector4d controller;
   	/// array of recorded skeletons
    std::vector< std::vector<Eigen::Matrix4d> > frameset;
	/// RBF interpolator handler
	RbfSkelInterp rbfint;
	/// compact support radius
	double csr;

	/**
	 * function for debugging purposes
	 */ 
	void printMarkers();

    public:
    
    /**
     * creates a new animation
     */
    Animation() {
		csr = 1;
	}
	
	/**
	 * clears all animation data
	 */
    ~Animation() {
        clearAnimation();
	}

    void setContext();
    Skeleton & getSkeleton();
    void setSkeleton(Skeleton & s);
    Eigen::Vector4d &getController();
    void setController(unsigned int fid);
    void setController(Eigen::Vector4d const & c);
    const std::vector< std::vector<Eigen::Matrix4d> > &getFrameSet();
    std::vector<Eigen::Vector4d> &getMarkers();
    void addNewKeyFrame(unsigned int fid, Skeleton &sk);
    void addNewKeyFrame(Eigen::Vector4d const & m, Skeleton &sk);

    /**
     * \brief clear all recorded poses
     */
    void clearAnimation() {
        markers.clear();
        frameset.clear();
    }
    /**
     * set compact support radius
     * \param csr new compact support radius
     */
    void setCsr(double r) {
		csr = r;
	}
	void animate();
};

#endif
