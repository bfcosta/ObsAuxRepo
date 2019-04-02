#ifndef KFSELECTOR
#define KFSELECTOR

#include <kfinterpolator.h>
#include <kf2dinterpolator.h>

class USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm
	 * @brief A base class for handling uniform sampling key frame selection 
	 */
	protected:
		/// desired number of keyframes
		unsigned int nkf;
		/// total of frames
		unsigned int tf;
		/// elected keyframes
		Eigen::VectorXi keyframes;
		/// elapsed time for selecting keyframes
		double eltime;
		/// Mocap interpolator handler
		MocapInterpolator * ih;

		USBaseSelector() {}
		/**
		 * @brief creates a simple uniform sampling algorithm for keyframe selection
		 * @param mi Mocap interpolator
		 * @param p percentage of keyframes
		 */		
		USBaseSelector(MocapInterpolator * mi, float p) : ih(mi) {
			tf = mi->getNumberOfFrames();
			nkf = tf * p;
			if (nkf < 3) nkf = 3;
		}
		/**
		 * @brief set interpolator strategy and percentage of keyframes
		 * @param mi Mocap interpolator
		 * @param p percentage of keyframes
		 */				
		void setInterpolator(MocapInterpolator * mi, float p) {
			ih = mi;
			tf = mi->getNumberOfFrames();
			nkf = tf * p;
			if (nkf < 3) nkf = 3;
		}
		virtual ~USBaseSelector() {}

	public:
		
		Eigen::VectorXi getKeyframeVector();
		virtual Eigen::VectorXi getKeyframes(unsigned int nkf);
		/**
		 * @brief get markers position for computed keyframes
		 * @param fids frame ids vector
		 * @return rowwise matrix for markers position vectors
		 */
		Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids) {
			return ih->getMarkers(fids);
		}
		/**
		 * @brief get total of keyframes
		 * @return number of keyframes
		 */
		unsigned int getNumberOfKeyframes() {
			return nkf;
		}
		/**
		 * @brief get frame distance matrix
		 * @return set frame distance matrix
		 */
		Eigen::MatrixXd getFrameDistanceMatrix() {
			return ih->getFrameDistanceMatrix();
		}
		/**
		 * @brief set frame distance matrix
		 * @param dmat frame distance matrix
		 */
		void setFrameDistanceMatrix(Eigen::MatrixXd & dmat) {
			ih->setFrameDistanceMatrix(dmat);
		}
		/**
		 * @brief get frame projection matrix
		 * @return set frame projection matrix
		 */
		Eigen::MatrixXd getFramesProjection() {
			return ih->getRdframes();
		}
		/**
		 * @brief set frame projection matrix
		 * @param rdmat frame projection matrix
		 */
		void setFramesProjection(Eigen::MatrixXd & rdmat) {
			ih->setRdframes(rdmat);
		}
		/**
		 * @brief get interpolated pose
		 * @param fid pose frame id
		 * @return skeleton pose
		 */
		Skeleton getInterpolatedPose(unsigned int fid) {
			return ih->getIFramePose(fid);
		}
		void printStats();
};

class PBBaseSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for position-based algorithm
	 * @brief A base class for handling position-based key frame selection 
	 */
	protected:
		unsigned int getMinErrorFrameId(const Eigen::VectorXd & ev, std::vector<bool> & iskf);
		PBBaseSelector() : USBaseSelector() {}
		/**
		 * @brief creates a simple position-based algorithm for keyframe extraction
		 * @param mocap interpolator
		 * @param p percentage of keyframes
		 */		
		PBBaseSelector(MocapInterpolator * mi, float p) : USBaseSelector(mi, p) {}
		virtual ~PBBaseSelector() {}

	public:
		virtual Eigen::VectorXi getKeyframes(unsigned int nkf);
		Eigen::VectorXi getKeyframesByNeigborhood(unsigned int nd = 5);
};

class SCSBaseSelector : public USBaseSelector  {
	/**
	 * @class keyframe extractor for SCS algorithm
	 * @brief A base class for handling SCS algorithm
	 */
	protected:
		int getMostDissimilarFrame(Eigen::VectorXi & kfs, int start, int len);
		/**
		 * @brief creates a simple SCS algorithm for keyframe extraction
		 * @param in mocap session handler
		 * @param mi mocap interpolator
		 * @param p percentage of keyframes
		 */
		SCSBaseSelector() : USBaseSelector() {}
		SCSBaseSelector(MocapInterpolator * mi, float p) : USBaseSelector(mi, p) {}
		virtual ~SCSBaseSelector() {}

	public:
		virtual Eigen::VectorXi getKeyframes(unsigned int nkf);
		Eigen::VectorXi getKeyframesByDeltaError(float ep = 0.9f);
};

class ROLSBaseSelector : public USBaseSelector  {
	/**
	 * @class keyframe extractor for ROLS algorithm
	 * @brief A base class for handling Regularized Orthogonal Least 
	 * Squares (ROLS) algorithm
	 */
	protected:
        /// RBF phi matrix
        Eigen::MatrixXd phimat;
        /// projection of subset data
        Eigen::MatrixXd ssproj;
        /// array of elected sample indices
        Eigen::VectorXi sample;

		/**
		 * @brief creates a simple ROLS algorithm for keyframe extraction
		 * @param in mocap session handler
		 * @param mi mocap interpolator
		 * @param p percentage of keyframes
		 */
		ROLSBaseSelector() : USBaseSelector() {}
		//ROLSBaseSelector(MocapInterpolator * mi, float p) :
		ROLSBaseSelector(RBFInterpolator * mi, float p) :
			USBaseSelector(mi, p) {}
		virtual ~ROLSBaseSelector() {}
		Eigen::VectorXi getSubsetData(double pc = 0.5);
		double phi(double d2);

	public:
		virtual Eigen::VectorXi getKeyframes(unsigned int nkf);
		Eigen::VectorXi chooseControlPoints(unsigned int nkf, double reg = 0.00001);
};

class ROLSRBFPBaseSelector : public ROLSBaseSelector  {
	/**
	 * @class keyframe extractor for ROLS algorithm with RBFP projection
	 * @brief A base class for handling Regularized Orthogonal Least 
	 * Squares (ROLS) algorithm with RBFP projection
	 */
	protected:

		/**
		 * @brief creates a simple ROLS algorithm  with RBFP projection for keyframe extraction
		 * @param in mocap session handler
		 * @param mi mocap interpolator
		 * @param p percentage of keyframes
		 */
		ROLSRBFPBaseSelector() : ROLSBaseSelector() {}
		ROLSRBFPBaseSelector(RBFInterpolator * mi, float p) :
			ROLSBaseSelector(mi, p) {}
		virtual ~ROLSRBFPBaseSelector() {}
		double projPhi(double d2);

	public:
		virtual Eigen::VectorXi getKeyframes(unsigned int nkf);
		Eigen::MatrixXd getRBF2DProjection(Eigen::VectorXi & kfs);
};

class FwdBkwdBaseSelector : public ROLSRBFPBaseSelector  {
	/**
	 * @class keyframe extractor algorithm with RBFP projection
	 * @brief A base class for handling a Forward and Backward selection
	 * algorithm with RBFP projection
	 */
	protected:

		/**
		 * @brief creates a simple FwdBkwd algorithm with RBFP projection for keyframe extraction
		 * @param in mocap session handler
		 * @param mi mocap interpolator
		 * @param p percentage of keyframes
		 */
		FwdBkwdBaseSelector() : ROLSRBFPBaseSelector() {}
		FwdBkwdBaseSelector(RBF2DInterpolator * mi, float p) :
			ROLSRBFPBaseSelector(mi, p) {}
		virtual ~FwdBkwdBaseSelector() {}

	public:
		virtual Eigen::VectorXi getKeyframes(unsigned int nkf);
};

/***********************************************************************/

class USEucLinSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * linear interpolation and positional DOFs
	 * @brief A class for handling uniform sampling key frame selection 
	 * with linear interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// linear mocap interpolator
		LinearInterpolator mci;
		
	public:
		USEucLinSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucSlerpSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spherical linear interpolation and positional DOFs
	 * @brief A class for handling uniform sampling key frame selection 
	 * with spherical linear interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// spherical linear mocap interpolator
		SphericalLinearInterpolator mci;
		
	public:
		USEucSlerpSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};


class USEucHerSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * hermitian interpolation and positional DOFs
	 * @brief A class for handling uniform sampling key frame selection 
	 * with hermitian interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// hermitian mocap interpolator
		HermiteInterpolator mci;
		
	public:
		USEucHerSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);	
			}
};

class USEucRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling uniform sampling key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF mocap interpolator
		RBFInterpolator mci;
		
	public:
		USEucRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);	
			}
};

class USEucForceRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by Force approach
		ForceRBFInterpolator mci;
		
	public:
		USEucForceRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucMdsRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by MDS approach
		MDSRBFInterpolator mci;
		
	public:
		USEucMdsRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucPcaRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by PCA approach
		PCARBFInterpolator mci;
		
	public:
		USEucPcaRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucLleRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LLE approach
		LLERBFInterpolator mci;
		
	public:
		USEucLleRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucLampRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		LAMPRBFInterpolator mci;
		
	public:
		USEucLampRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucTsneRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by t-SNE approach
		TsneRBFInterpolator mci;
		
	public:
		USEucTsneRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class USEucIsomapRbfSelector : public USBaseSelector {
	/**
	 * @class keyframe selector for uniform sampling algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by isomap approach
		IsomapRBFInterpolator mci;
		
	public:
		USEucIsomapRbfSelector(MocapReader * mr, float p = 0.1f) : 
			USBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class PBEucLinSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * linear interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with linear interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// linear mocap interpolator
		LinearInterpolator mci;
		
	public:
		PBEucLinSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class PBEucSlerpSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spherical linear interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spherical linear interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// spherical linear mocap interpolator
		SphericalLinearInterpolator mci;
		
	public:
		PBEucSlerpSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class PBEucHerSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * hermitian interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with hermitian interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// hermitian mocap interpolator
		HermiteInterpolator mci;
		
	public:
		PBEucHerSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class PBEucRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator
		RBFInterpolator mci;
		
	public:
		PBEucRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class PBEucForceRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by Force approach
		ForceRBFInterpolator mci;
		
	public:
		PBEucForceRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class PBEucMdsRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by MDS approach
		MDSRBFInterpolator mci;
		
	public:
		PBEucMdsRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class PBEucPcaRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by PCA approach
		PCARBFInterpolator mci;
		
	public:
		PBEucPcaRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class PBEucLleRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by LLE approach
		LLERBFInterpolator mci;
		
	public:
		PBEucLleRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class PBEucLampRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by LAMP approach
		LAMPRBFInterpolator mci;
		
	public:
		PBEucLampRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class PBEucTsneRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by t-SNE approach
		TsneRBFInterpolator mci;
		
	public:
		PBEucTsneRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class PBEucIsomapRbfSelector : public PBBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * RBF/SKF and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with RBF/SKF scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with makers spread by isomap approach
		IsomapRBFInterpolator mci;
		
	public:
		PBEucIsomapRbfSelector(MocapReader * mr, float p = 0.1f) : 
			PBBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucLinSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * linear interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with linear interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// linear mocap interpolator
		LinearInterpolator mci;
		
	public:
		SCSEucLinSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class SCSEucSlerpSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spherical linear interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spherical linear interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// spherical linear mocap interpolator
		SphericalLinearInterpolator mci;
		
	public:
		SCSEucSlerpSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class SCSEucHerSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * hermitian interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with hermitian interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// hermitian mocap interpolator
		HermiteInterpolator mci;
		
	public:
		SCSEucHerSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class SCSEucRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF mocap interpolator
		RBFInterpolator mci;
		
	public:
		SCSEucRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class SCSEucForceRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by Force approach
		ForceRBFInterpolator mci;
		
	public:
		SCSEucForceRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucMdsRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by MDS approach
		MDSRBFInterpolator mci;
		
	public:
		SCSEucMdsRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucPcaRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by PCA approach
		PCARBFInterpolator mci;
		
	public:
		SCSEucPcaRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucLleRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LLE approach
		LLERBFInterpolator mci;
		
	public:
		SCSEucLleRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucLampRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		LAMPRBFInterpolator mci;
		
	public:
		SCSEucLampRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucTsneRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by t-SNE approach
		TsneRBFInterpolator mci;
		
	public:
		SCSEucTsneRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class SCSEucIsomapRbfSelector : public SCSBaseSelector {
	/**
	 * @class keyframe selector for position based algorithm with
	 * spatial keyframe interpolation and positional DOFs
	 * @brief A class for handling position based key frame selection 
	 * with spatial keyframe interpolation scheme and euclidean distance for DOFs
	 */
	protected:
		/// skeleton analizer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by isomap approach
		IsomapRBFInterpolator mci;
		
	public:
		SCSEucIsomapRbfSelector(MocapReader * mr, float p = 0.1f) : 
			SCSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, single dimension
	 * markers and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, single dimension markers and
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with 1D markers
		RBFInterpolator mci;

	public:
		RolsEucRbfSelector(MocapReader * mr, float p = 0.1f) :
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class RolsRbfpEucRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, single dimension
	 * markers and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, single dimension, markers and
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with 1D markers
		RBFInterpolator mci;

	public:
		RolsRbfpEucRbfSelector(MocapReader * mr, float p = 0.1f) :
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class RolsEucForceRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, Force dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, Force dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by Force approach
		ForceRBFInterpolator mci;
		
	public:
		RolsEucForceRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucForceRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, Force dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, Force dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by Force approach
		ForceRBFInterpolator mci;
		
	public:
		RolsRbfpEucForceRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucForceRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward 
	 * algorithm with spatial keyframe interpolation, Force dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, Force dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by Force approach
		ForceRBFInterpolator mci;
		
	public:
		FwdBkwdEucForceRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucMdsRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, MDS dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, MDS dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by MDS approach
		MDSRBFInterpolator mci;
		
	public:
		RolsEucMdsRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucMdsRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, MDS dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, MDS dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by MDS approach
		MDSRBFInterpolator mci;
		
	public:
		RolsRbfpEucMdsRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucMdsRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward
	 * algorithm with spatial keyframe interpolation, MDS dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, MDS dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by MDS approach
		MDSRBFInterpolator mci;
		
	public:
		FwdBkwdEucMdsRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucPcaRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, PCA dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, MDS dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by PCA approach
		PCARBFInterpolator mci;
		
	public:
		RolsEucPcaRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucPcaRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, PCA dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, MDS dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by PCA approach
		PCARBFInterpolator mci;
		
	public:
		RolsRbfpEucPcaRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucPcaRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward
	 * algorithm with spatial keyframe interpolation, PCA dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, MDS dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by PCA approach
		PCARBFInterpolator mci;
		
	public:
		FwdBkwdEucPcaRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucLleRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, LLE dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LLE approach
		LLERBFInterpolator mci;
		
	public:
		RolsEucLleRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucLleRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, LLE dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LLE approach
		LLERBFInterpolator mci;
		
	public:
		RolsRbfpEucLleRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucLleRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward
	 * algorithm with spatial keyframe interpolation, LLE dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LLE approach
		LLERBFInterpolator mci;
		
	public:
		FwdBkwdEucLleRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucTsneRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, t-SNE dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, t-SNE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by t-SNE approach
		TsneRBFInterpolator mci;
		
	public:
		RolsEucTsneRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucTsneRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, t-SNE dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, t-SNE dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by t-SNE approach
		TsneRBFInterpolator mci;
		
	public:
		RolsRbfpEucTsneRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucTsneRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward
	 * algorithm with spatial keyframe interpolation, t-SNE dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, t-SNE dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by t-SNE approach
		TsneRBFInterpolator mci;
		
	public:
		FwdBkwdEucTsneRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucIsomapRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, isomap dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, isomap dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by isomap approach
		IsomapRBFInterpolator mci;
		
	public:
		RolsEucIsomapRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucIsomapRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, isomap dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, isomap dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by isomap approach
		IsomapRBFInterpolator mci;
		
	public:
		RolsRbfpEucIsomapRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucIsomapRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward
	 * algorithm with spatial keyframe interpolation, isomap dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, isomap dimension reduction and 
	 * euclidean distance for DOFs. 2D projection is done by RBFP.
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by isomap approach
		IsomapRBFInterpolator mci;
		
	public:
		FwdBkwdEucIsomapRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucLampRbfSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, LAMP dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		LAMPRBFInterpolator mci;
		
	public:
		RolsEucLampRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsRbfpEucLampRbfSelector : public ROLSRBFPBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, LAMP dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		LAMPRBFInterpolator mci;
		
	public:
		RolsRbfpEucLampRbfSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSRBFPBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class FwdBkwdEucLampRbfSelector : public FwdBkwdBaseSelector {
	/**
	 * @class keyframe selector for forward and backward
	 * algorithm with spatial keyframe interpolation, LAMP dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling forward and backward key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		LAMPRBFInterpolator mci;
		
	public:
		FwdBkwdEucLampRbfSelector(MocapReader * mr, float p = 0.1f) : 
			FwdBkwdBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
		/**
		 * @brief get the merged 2D projection with another mocap session
		 * @param sib sibling mocap session
		 * @return 2D projection
		 */
		Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) {
			return mci.getMerged2DProjection(sib);
		}
};

class RolsEucLinSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, LAMP dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		LinearInterpolator mci;
		
	public:
		RolsEucLinSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

class RolsEucHerSelector : public ROLSBaseSelector {
	/**
	 * @class keyframe selector for regularized orthogonal least squares
	 * algorithm with spatial keyframe interpolation, LAMP dimension 
	 * reduction approach and positional DOFs
	 * @brief A class for handling ROLS key frame selection with spatial
	 * keyframe interpolation scheme, LLE dimension reduction and 
	 * euclidean distance for DOFs
	 */
	protected:
		/// skeleton analyzer for euclidean distance
		GeodesicSkeletonAnalizer esa;
		/// RBF/SKF mocap interpolator with markers spread by LAMP approach
		HermiteInterpolator mci;
		
	public:
		RolsEucHerSelector(MocapReader * mr, float p = 0.1f) : 
			ROLSBaseSelector(), esa(), mci(mr, &esa) {
				setInterpolator(&mci, p);
			}
};

#endif
