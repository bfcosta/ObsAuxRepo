#ifndef KFINTERPOLATOR
#define KFINTERPOLATOR

#include <skanalizer.h>
#include <mocapreader.h>
#include <animation.h>

class MocapInterpolator {
	/**
	 * @class Interface for handling interpolation algorithms
	 * @brief A base class to handle different interpolation strategies for Mocap rebuild.
	 */
	protected:
		/// skeleton distance helper
		SkeletonAnalizer * skdist;
		/// mocap session handler
		MocapReader * mc;
		/// frame error vector
		Eigen::VectorXd evec;
		/// frame distance matrix
		Eigen::MatrixXd fdist;
		/// frame pose in reduced dimension (for 2D)
		Eigen::MatrixXd rdframes;
		
	public:
		MocapInterpolator(MocapReader * mr, SkeletonAnalizer * sa) {
			mc = mr;
			skdist = sa;
		}
		virtual ~MocapInterpolator() {}
		virtual Eigen::VectorXd rebuildMocapData(Eigen::VectorXi & kfs, int start, int len) = 0;
		virtual double getFrameDistance(unsigned int pf, unsigned int ff, unsigned int nf) = 0;
		virtual double getFrameDistance(Eigen::VectorXi & kfs, int fid) = 0;
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids) = 0;
		virtual void setRdframes() = 0;
		virtual Eigen::VectorXi getINeighbourFrames(Eigen::VectorXi const & kfs, unsigned int fid) = 0;
		virtual Eigen::VectorXd getProjectionErrorVector(Eigen::VectorXi & kfs, int start, int len) = 0;
		virtual double getProjectedFrameDistance(Eigen::VectorXi & kfs, int fid) = 0;
		virtual Skeleton getIFramePose(unsigned int fid) = 0;
		
		/**
		 * @brief get average reconstruction error
		 * @return average error
		 */
		double getAvgError() {
			double r = mc->getNumberOfDOFs() * mc->getNumberOfFrames();
			return (evec.sum() / r);	
		}
		/**
		 * @brief return previous error calculated
		 * @return error vector
		 */
		Eigen::VectorXd getErrorVector() {
			return evec;
		}
		/**
		 * @brief return total of frames from mocap session
		 * @return total of frames
		 */
		unsigned int getNumberOfFrames() {
			return mc->getNumberOfFrames();
		}
		/**
		 * @brief set frame projection matrix
		 * @param frame projection matrix
		 */
		void setRdframes(Eigen::MatrixXd & rdmat) {
			rdframes = rdmat;
		}
		/**
		 * @brief get current frame projection matrix
		 * @return frame projection matrix
		 */
		Eigen::MatrixXd getRdframes() {
			return rdframes;
		}
		
		Eigen::MatrixXd getDistanceMatrix(Eigen::VectorXi const & fref);
		Eigen::MatrixXd getFrameDistanceMatrix();
		void setFrameDistanceMatrix(Eigen::MatrixXd & dmat);
};

class LinearInterpolator: public MocapInterpolator {
	/**
	 * @class Class for handling linear interpolation
	 * @brief A class for handling linear interpolation of Mocap data
	 */
	protected:
		/// mocap data in matrix format, rotation given in euler angles
		Eigen::MatrixXd dofs;
		/// rebuilt data in matrix format, rotation given in euler angles
		Eigen::MatrixXd idofs;

	public:
		/**
		 * @brief buids a mocap linear interpolator
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		LinearInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : MocapInterpolator(mr, sa) {
			dofs = mr->getEulerAnglesDOFMatrix();
			idofs = dofs;
			evec = Eigen::ArrayXd::Zero(dofs.rows());
		}
		
		virtual ~LinearInterpolator() {}
		virtual Eigen::VectorXd rebuildMocapData(Eigen::VectorXi & kfs, int start, int len);
		virtual double getFrameDistance(unsigned int pf, unsigned int ff, unsigned int nf);
		virtual double getFrameDistance(Eigen::VectorXi & kfs, int fid);
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::VectorXi getINeighbourFrames(Eigen::VectorXi const & kfs, unsigned int fid);
		virtual Skeleton getIFramePose(unsigned int fid);
		/**
		 * @brief hidden method for handling 2D markers with RBF interpolator.
		 * Does not do anything in this class.
		 */
		virtual void setRdframes() {}
		/**
		 * @brief get the projected frames error vector for given keyframes and interval specified
		 * @param kfs current keyframes
		 * @param start frame id to start rebuild
		 * @param len number of frames to rebuild
		 * @return projected frame error vector
		 */
		virtual Eigen::VectorXd getProjectionErrorVector(Eigen::VectorXi & kfs, int start, int len) {
			return rebuildMocapData(kfs, start, len);
		}
		/**
		 * @brief get the projected frame distance between a marker given by frame id
		 * and the current interpolation given by a key frame set. The key frame
		 * set is given as an ordered vector.
		 * @param kfs key frames id set as an ordered vector
		 * @param fid projected frame id to be evaluated
		 * @return pose distance between projected frame and interpolation reconstruction
		 */
		virtual double getProjectedFrameDistance(Eigen::VectorXi & kfs, int fid) {
			return getFrameDistance(kfs, fid);
		}
};

class SphericalLinearInterpolator: public LinearInterpolator {
	/**
	 * @class Class for handling spherical linear interpolation
	 * @brief A class for handling spherical linear interpolation of Mocap data
	 */
	public:
		/**
		 * @brief buids a mocap spherical linear interpolator
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		SphericalLinearInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : LinearInterpolator(mr, sa) {
		}
		
		virtual ~SphericalLinearInterpolator() {}
		virtual double getFrameDistance(unsigned int pf, unsigned int ff, unsigned int nf);
};

class HermiteInterpolator: public LinearInterpolator {
	/**
	 * @class Class for handling cubic hermitian interpolation
	 * @brief A class for handling cubic hermitian interpolation of Mocap data
	 */
	protected:
		/// first derivative of mocap data in matrix format, rotation given in euler angles
		Eigen::MatrixXd d1dofs;
		
		void setDerivative();
		
	public:
		/**
		 * @brief buids a mocap hermitian interpolator
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		HermiteInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : LinearInterpolator(mr, sa) {
			setDerivative();
		}
		virtual ~HermiteInterpolator() {}
		virtual double getFrameDistance(unsigned int pf, unsigned int ff, unsigned int nf);
		//virtual double getFrameDistance(Eigen::VectorXi & kfs, unsigned int fid);
};

class RBFInterpolator: public MocapInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data
	 */
	protected:
		/// skeleton pose set for equivalent key frame
		std::vector<Skeleton> poseset;
		/// skf animation handler
		Animation animHandler;
		/// controller position for animation
		Eigen::Vector4d controller;		
		/// spatial markers for keyframed pose
		std::vector<Eigen::Vector4d> markers;
		
		double getAnimationError(Eigen::MatrixXd const & marker, Eigen::VectorXi const & kfs, int ff);
		double getAnimationError(Eigen::MatrixXd const & marker, int ff);
		double getCSRadius(Eigen::VectorXi & kfs);
		void addTimeCoordinate(Eigen::MatrixXd & proj);
		Eigen::MatrixXd getSelectedRdframes(Eigen::VectorXi const & fids);
		
	public:
		/**
		 * @brief buids a mocap RBF interpolator
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		RBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : MocapInterpolator(mr, sa), animHandler() {
			evec = Eigen::ArrayXd::Zero(mr->getNumberOfFrames());
		}
		
		virtual ~RBFInterpolator() {}
		virtual Eigen::VectorXd rebuildMocapData(Eigen::VectorXi & kfs, int start, int len);
		virtual double getFrameDistance(unsigned int pf, unsigned int ff, unsigned int nf);
		virtual double getFrameDistance(Eigen::VectorXi & kfs, int fid);
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::VectorXi getINeighbourFrames(Eigen::VectorXi const & kfs, unsigned int fid);
		virtual Skeleton getIFramePose(unsigned int fid);
		virtual void setRdframes();
		void rescaleMarkers(Eigen::MatrixXd & proj);
		/**
		 * @brief get the projected frames error vector for given keyframes and interval specified
		 * @param kfs current keyframes
		 * @param start frame id to start rebuild
		 * @param len number of frames to rebuild
		 * @return projected frame error vector
		 */
		virtual Eigen::VectorXd getProjectionErrorVector(Eigen::VectorXi & kfs, int start, int len) {
			return rebuildMocapData(kfs, start, len);
		}
		/**
		 * @brief get the projected frame distance between a marker given by frame id
		 * and the current interpolation given by a key frame set. The key frame
		 * set is given as an ordered vector.
		 * @param kfs key frames id set as an ordered vector
		 * @param fid projected frame id to be evaluated
		 * @return pose distance between projected frame and interpolation reconstruction
		 */
		virtual double getProjectedFrameDistance(Eigen::VectorXi & kfs, int fid) {
			return getFrameDistance(kfs, fid);
		}
};

#endif
