#ifndef BVHREADER_H
#define BVHREADER_H

#include <vector>
#include <mocapreader.h>

/**
 * \class BVHReader
 *
 * \brief A reader for BVH file format.
 *
 * Implements a reader to BVH file format. This format has two sections,
 * both inside the same file: a skeleton description on its resting pose
 * and data for input channels to change the pose. On its
 * first section called HIERARCHY, most important tags are OFFSET and 
 * CHANNELS. The first describes the bone length and the second is an 
 * interface for merging data provided at the second section called MOTION. 
 * Data is provided on depth-first search at MOTION section. Details are
 * omitted here. 
 * 
 * \see http://research.cs.wisc.edu/graphics/Courses/cs-838-1999/Jeff/BVH.html
 */

class BVHReader: public MocapReader {
	protected:
	
	/**
	 * A structure to read OFFSET tag on BVH file
	 */
	struct triplef {
		float x,y,z;
	};
	
	/**
	 * A structure to read CHANNELS tag on BVH file.
	 */
	struct mask {
		/// array of rotation/translation order, X, Y or Z.
		char dim[3];
		/// flag to describe a rotational/translational channel
		bool rot;
		/// offset joint reference
		unsigned int ref;
	};

	/// bone name at file
	std::vector<std::string> name;
	/// OFFSET value at file
	std::vector<struct triplef> joint;
	/// sequence of channels found on BVH file
	std::vector<struct mask> channel;
	/// parent joint
	std::vector<unsigned int> parent;
	/// number of frames
	unsigned int nframe;
	/// frame time interval
	float tframe;
	/// bone maximum length
	float mblen;
	/// array of channel sequence values describing all poses
	std::vector< std::vector<float> > motion;
	void readBvh(const std::string & fname);

	public:
	BVHReader();
	virtual void readFile(const std::string & fname);
	virtual Skeleton getSkeleton();
	virtual unsigned int getNumberOfDOFs();
	virtual unsigned int getNumberOfFrames();
	virtual Skeleton getFrameSkeleton(unsigned int i);
	virtual Skeleton getFrameSkeleton(Eigen::VectorXd dofs, bool euler = true);
	virtual unsigned int getFrameInterval();
	virtual void print();
	virtual Eigen::MatrixXd getEulerAnglesDOFMatrix();
	virtual Eigen::MatrixXd getQuaternionDOFMatrix();
	virtual void printPose(Skeleton & sk, std::ostream & ost);
	virtual Eigen::VectorXd getEulerAnglesDOFPose(Skeleton & sk);
};

#endif
