#ifndef MOCAPREADER_H
#define MOCAPREADER_H

#include <skeleton.h>
#include <string>
#include <ostream>

/**
 * \class MocapReader
 *
 * \brief Interface class for dealing with many animation file formats
 *
 * Interface for reading animation file formats. It should be able to
 * read an animation file and a skeleton in its resting pose or in a
 * specific frame pose, if the file format has been implemented.
 * Details are on its subclasses.
 * 
 * \see BVHReader
 * \see ASFReader
 */

class MocapReader {
	protected:
		/// thickness ratio between bone length and cylinder radius
		static const float lcratio = 25;

		/**
		 * @brief break rotation matrix in euler angles where rotation order is ZXY
		 * @param rot rotation matrix
		 * @return vector with euler angles
		 */
		Eigen::Vector3d getZXYRotation(Eigen::Matrix4d & rot) {
			Eigen::Vector3d ret;
			double cos0;

			// rotation order is ZXY
			// get x rotation
			ret(0) = std::asin(rot(2,1));
			if (ret(0) > M_PI / 2)
				ret(0) = M_PI - ret(0);
			cos0 = std::cos(ret(0));
			if ((rot(2,1) == 1) || (rot(2,1) == -1)) {
				// get y rotation, z = 0
				ret(1) = std::atan2(rot(1,0), rot(0,0));
			} else {
				// get y and z rotation
				ret(2) = std::atan2(-1*rot(0,1) / cos0 , rot(1,1) / cos0);
				ret(1) = std::atan2(-1*rot(2,0) / cos0 , rot(2,2) / cos0);
			}
			ret *= 180/M_PI;
			return ret;
		}

		/**
		 * @brief break rotation matrix in euler angles where rotation order is ZYX
		 * @param rot rotation matrix
		 * @return vector with euler angles
		 */
		Eigen::Vector3d getZYXRotation(Eigen::Matrix4d & rot) {
			Eigen::Vector3d ret;
			double cos0;

			// rotation order is ZYX
			// get y rotation
			ret(1) = std::asin(-1*rot(2,0));
			if (ret(1) > M_PI / 2)
				ret(1) = M_PI - ret(1);
			cos0 = std::cos(ret(1));
			if ((rot(2,0) == 1) || (rot(2,0) == -1)) {
				// get x rotation, z = 0
				ret(0) = std::atan2(rot(0,1), rot(0,2));
			} else {
				// get x and z rotation
				ret(2) = std::atan2(rot(1,0) / cos0 , rot(0,0) / cos0);
				ret(0) = std::atan2(rot(2,1) / cos0 , rot(2,2) / cos0);
			}
			ret *= 180/M_PI;
			return ret;
		}

	
	public:
		/**
		 * clears all data
		 */
		virtual ~MocapReader() {}
		/**
		 * \param fname animation file
		 * 
		 * reads file and loads animation data
		 */
		virtual void readFile(const std::string & fname) = 0;
		/**
		 * \return skeleton on resting pose
		 * 
		 * get the skeleton resting pose
		 */
		virtual Skeleton getSkeleton() = 0;
		/**
		 * \return number of DOFs in animated skeleton
		 * 
		 * get number of DOFs in animated skeleton
		 */
		virtual unsigned int getNumberOfDOFs() = 0;
		/**
		 * \return number of frames in animation
		 * 
		 * get number of frames in animation
		 */
		virtual unsigned int getNumberOfFrames() = 0;
		/**
		 * \param i frame id counter
		 * \return a skeleton pose
		 * 
		 * returns a skeleton at a specific frame id during animation
		 */
		virtual Skeleton getFrameSkeleton(unsigned int i) = 0;
		/**
		 * \param i dof configuration in vector format
		 * \param euler format of mocap matrix. True is euler angles and false 
		 * is quaternion.
		 * \return a skeleton pose
		 * 
		 * returns a skeleton at a specific pose
		 */
		virtual Skeleton getFrameSkeleton(Eigen::VectorXd dofs, bool euler = true) = 0;
		/**
		 * \return frame time interval
		 * 
		 * get frame time interval in milisseconds
		 */
		virtual unsigned int getFrameInterval() = 0;
		/**
		 * dump structure to stdout for debugging purposes
		 */
		virtual void print() = 0;
		/**
		 * \return mocap DOF matrix with loaded mocap data
		 * 
		 * dump all frames as one matrix. Rows are frames and columns are 
		 * skeleton DOFs. Joint data are described as euler angles.
		 */
		virtual Eigen::MatrixXd getEulerAnglesDOFMatrix() = 0;
		/**
		 * \return mocap DOF matrix with loaded mocap data
		 * 
		 * dump all frames as one matrix. Rows are frames and columns are 
		 * skeleton DOFs. Joint data are described as quaternions.
		 */
		virtual Eigen::MatrixXd getQuaternionDOFMatrix() = 0;
		/**
		 * \param sk skeleton pose
		 * \param ost output stream to print
		 * 
		 * print skeleton pose as a current format frame in file
		 */
		virtual void printPose(Skeleton & sk, std::ostream & ost) = 0;
		/**
		 * \return mocap DOF vector with skeleton pose data
		 * 
		 * dump skeleton pose as a DOF vector. Joint data are described as euler angles.
		 */
		virtual Eigen::VectorXd getEulerAnglesDOFPose(Skeleton & sk) = 0;
};

#endif
