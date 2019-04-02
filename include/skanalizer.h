#ifndef SKANALIZER_H
#define SKANALIZER_H

#include <Eigen/Dense>
#include <skeleton.h>

/**
 * @class
 * @brief
 */
class SkeletonAnalizer {
	protected:
		/// weight function flag
		unsigned int wfunc;
		std::vector<double> getWeights(Skeleton & s1, unsigned int option = 0);

	public:
		/**
		 * @brief
		 */
		SkeletonAnalizer() {}
		/**
		 * @brief
		 */
		virtual ~SkeletonAnalizer() {}
		/**
		 * @brief
		 */
		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2) = 0;
		/**
		 * @brief
		 */
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1) = 0;
};

/**
 * @class
 * @brief
 */
class EuclideanSkeletonAnalizer: public SkeletonAnalizer {

	public:
		/**
		 * @brief builds a skeleton distance helper for euclidean distance metric.
		 * @param w weight function flag. 0 = constant, 1 = increasing,
		 * 2 = decreasing.
		 */
		EuclideanSkeletonAnalizer(unsigned int w = 0) : SkeletonAnalizer() {
			wfunc = w;
		}
		
		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
		double getMatRotDist(Eigen::Matrix4d const & m1, Eigen::Matrix4d const & m2);
		double getMatEuclDist(Eigen::Matrix4d const & m1, Eigen::Matrix4d const & m2);
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1);
};

/**
 * @class
 * @brief
 */
class TrInvEucSkeletonAnalizer: public EuclideanSkeletonAnalizer {

	public:
		/**
		 * @brief builds a skeleton distance helper for translation invariant
		 * euclidean distance metric.
		 * @param w weight function flag. 0 = constant, 1 = increasing,
		 * 2 = decreasing.
		 */
		TrInvEucSkeletonAnalizer(unsigned int w = 0) : EuclideanSkeletonAnalizer(w) {
		}

		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1);
};

/**
 * @class
 * @brief
 */
class RootlessEucSkeletonAnalizer: public EuclideanSkeletonAnalizer {

	public:
		/**
		 * @brief builds a skeleton distance helper for euclidean distance
		 * metric without root joint.
		 * @param w weight function flag. 0 = constant, 1 = increasing,
		 * 2 = decreasing.
		 */
		RootlessEucSkeletonAnalizer(unsigned int w = 0) : EuclideanSkeletonAnalizer(w) {
		}

		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1);
};

/**
 * @class
 * @brief
 */
class QuaternionSkeletonAnalizer: public SkeletonAnalizer {
	protected:
		// weight function flag
		//unsigned int wfunc;
		//std::vector<double> getWeights(Skeleton & s1, unsigned int option = 0);
		virtual double getDistance(Eigen::Matrix4d const & m1, Eigen::Matrix4d const & m2);

	public:
		/**
		 * @brief buils a skeleton distance helper for quaternion distance metric.
		 * @param w weight function flag. 0 = constant, 1 = increasing,
		 * 2 = decreasing, 3 = bone length, 4 = hierarchical bone length.
		 */
		QuaternionSkeletonAnalizer(unsigned int w = 4) : SkeletonAnalizer() {
			wfunc = w;
		}

		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1);
};

/**
 * @class
 * @brief
 */
class RootlessQuatSkeletonAnalizer: public QuaternionSkeletonAnalizer {

	public:
		/**
		 * @brief buils a skeleton distance helper for quaternion distance metric without root joint
		 * @param w weight function flag. 0 = constant, 1 = increasing,
		 * 2 = decreasing, 3 = bone length, 4 = hierarchical bone length.
		 */
		RootlessQuatSkeletonAnalizer(unsigned int w = 4) : QuaternionSkeletonAnalizer(w) {
		}

		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1);
};


/**
 * @class
 * @brief
 */
class GeodesicSkeletonAnalizer: public QuaternionSkeletonAnalizer {
	protected:
		virtual double getDistance(Eigen::Matrix4d const & m1, Eigen::Matrix4d const & m2);

	public:
		/**
		 * @brief builds a skeleton distance helper for quaternion distance metric.
		 * @param w weight function flag. 0 = constant, 1 = increasing,
		 * 2 = decreasing, 3 = bone length, 4 = hierarchical bone length.
		 */
		GeodesicSkeletonAnalizer(unsigned int w = 5) : QuaternionSkeletonAnalizer(w) {}
		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
		virtual Eigen::VectorXd getMultidimensionalWeightedPoseVector(Skeleton & s1);
};

/**
 * @class
 * @brief
 */
class AreaSkeletonAnalizer: public EuclideanSkeletonAnalizer {
	protected:
		double getBoneArea(Eigen::Vector4d const & b00, Eigen::Vector4d const & b01, 
			Eigen::Vector4d const & b10, Eigen::Vector4d const & b11);

	public:
		/**
		 * @brief
		 */
		AreaSkeletonAnalizer(unsigned int w = 1) : EuclideanSkeletonAnalizer(w) {}
		virtual double getPoseDistance(Skeleton & s1, Skeleton & s2);
};

#endif
