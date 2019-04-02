#ifndef KF2DINTERPOLATOR
#define KF2DINTERPOLATOR

#include <kfinterpolator.h>

class RBF2DInterpolator: public RBFInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with 2D markers and
	 * aggregating common methods. Only sons of this class should be used.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */
	protected:
		double getProjectedFrameDistance(unsigned int pf, unsigned int ef, unsigned int nf);
		Eigen::MatrixXd getMergedDistanceMatrix(MocapReader * sib);

	public:
		/**
		 * @brief builds a mocap RBF interpolator with markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		RBF2DInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBFInterpolator(mr, sa) {
		}
		
		virtual ~RBF2DInterpolator() {}
		virtual Eigen::VectorXd getProjectionErrorVector(Eigen::VectorXi & kfs, int start, int len);
		virtual double getProjectedFrameDistance(Eigen::VectorXi & kfs, int fid);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib) = 0;
		Eigen::VectorXd getRebuildErrorVector(Eigen::MatrixXd const & proj, Animation * skfanim);
		Eigen::MatrixXd getGradientMatrix(Eigen::MatrixXd const & proj, Animation * skfanim, Eigen::VectorXd const & cerr, double alpha);
		Eigen::MatrixXd improveProjectionByGD(Eigen::MatrixXd const & proj, Eigen::VectorXi & keyframes, unsigned int iter = 200, double alpha = 0.5);

};

class ForceRBFInterpolator: public RBF2DInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with Force projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		Eigen::MatrixXd getDistanceMatrixProjection(Eigen::MatrixXd & ospdist, 
			unsigned int train = 50, double stepsize = 0.125);
		Eigen::MatrixXd getForceApproachProjection(Eigen::VectorXi const & fref, 
			unsigned int train = 50, double stepsize = 0.125);
		
	public:
		/**
		 * @brief builds a mocap RBF interpolator with Force strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		ForceRBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBF2DInterpolator(mr, sa) {
		}
		
		virtual ~ForceRBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};

class MDSRBFInterpolator: public RBF2DInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with MDS projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		Eigen::MatrixXd getDistanceMatrixProjection(Eigen::MatrixXd & adist);
		
	public:
		/**
		 * @brief buids a mocap RBF interpolator with MDS strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		MDSRBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBF2DInterpolator(mr, sa) {
		}
		
		virtual ~MDSRBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};

class IsomapRBFInterpolator: public RBF2DInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with Isomap projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		Eigen::MatrixXd getIsomapDistanceMatrix(Eigen::VectorXi const & fref, unsigned int nn = 14);
		Eigen::MatrixXd getIsomapProjection(Eigen::MatrixXd const & dist);
		Eigen::MatrixXd getMergedIsomapDistanceMatrix(MocapReader * sib, unsigned int nn = 14);
		void runFloydWarshallUpdate(Eigen::MatrixXd & fmat);
		bool sameMatrices(Eigen::MatrixXd & mat1, Eigen::MatrixXd & mat2);
		
	public:
		/**
		 * @brief buids a mocap RBF interpolator with Isomap strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		IsomapRBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBF2DInterpolator(mr, sa) {
		}
		
		virtual ~IsomapRBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};

class PCARBFInterpolator: public RBF2DInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with PCA projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		Eigen::MatrixXd getMergedRowwiseDataMatrix(MocapReader * sib);
		Eigen::MatrixXd getRowwiseDataMatrix(Eigen::VectorXi const & fids);
		Eigen::MatrixXd getPCAProjection(Eigen::MatrixXd const & datmat);
		
	public:
		/**
		 * @brief buids a mocap RBF interpolator with PCA strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		PCARBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBF2DInterpolator(mr, sa) {
		}
		
		virtual ~PCARBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};


class LLERBFInterpolator: public RBF2DInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with LLE projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		Eigen::VectorXi getKnnIndices(Eigen::VectorXd const & dvec, unsigned int knn, unsigned int curr);
		Eigen::MatrixXd getLLEWeightMatrix(Eigen::MatrixXd & dist, unsigned int knn);
		Eigen::MatrixXd getEmbeddingCoordinateMatrix(Eigen::MatrixXd const & wmat, unsigned int dim);
		Eigen::MatrixXd getDistanceMatrixProjection(Eigen::MatrixXd & dmat, unsigned int knn);
		
	public:
		/**
		 * @brief buids a mocap RBF interpolator with LLE strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		LLERBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBF2DInterpolator(mr, sa) {
		}
		
		virtual ~LLERBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};

class LAMPRBFInterpolator: public ForceRBFInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with LAMP projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		int getMostDistantFrameFromExtremes(Eigen::VectorXi const & fids);
		Eigen::VectorXi getForceSubset(Eigen::VectorXi const & ref, int ssetsize);
		Eigen::Vector4d getMedianProjection(unsigned int median, unsigned int nf);
		Eigen::VectorXd getLampProjection(Eigen::VectorXd const & mpose, 
		Eigen::MatrixXd const & yproj, Eigen::MatrixXd const & xproj);
			
	public:
		/**
		 * @brief buids a mocap RBF interpolator with LAMP strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		LAMPRBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : ForceRBFInterpolator(mr, sa) {
		}
		
		virtual ~LAMPRBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};

class TsneRBFInterpolator: public RBF2DInterpolator {
	/**
	 * @class Class for handling Spatial Keyframe RBF interpolation with t-SNE projection approach for 2D markers.
	 * @brief A class for handling Spatial Keyframe RBF interpolation of Mocap data.
	 */

	protected:
		Eigen::MatrixXd getRandomMatrix(unsigned int rows, unsigned int cols,
			double mu = 0, double sigma = 1);
		Eigen::VectorXd getProbabilitiesBySigma(Eigen::VectorXd const & dvec,
			double sigma, unsigned int pivot, double *h);
		Eigen::MatrixXd getPairwiseAffinities(Eigen::MatrixXd const & dmat,
			double perplexity = 30.0);
		Eigen::MatrixXd getLowDimensionalAffinities(Eigen::MatrixXd const & smat);
		Eigen::MatrixXd getKLGradient(Eigen::MatrixXd const & pmat,
			Eigen::MatrixXd const & qmat, Eigen::MatrixXd const & ymat);
		double getKLDivergence(Eigen::MatrixXd const & pmat,
			Eigen::MatrixXd const & qmat);
		Eigen::MatrixXd applyAndUpdateGain(Eigen::MatrixXd const & grad,
			Eigen::MatrixXd & gains, Eigen::MatrixXd const & deltaY, double ming = 0.01);
		void centralizeData(Eigen::MatrixXd & data);
		Eigen::MatrixXd getTsneDimensionReduction(Eigen::MatrixXd const & pij,
			Eigen::MatrixXd const & ymat, unsigned int train = 1000, double learn = 500,
			double momentum = 0.5);

	public:
		/**
		 * @brief builds a mocap RBF interpolator with t-SNE strategy for
		 * spreading markers at 2D space
		 * @param mr mocap session
		 * @param sa skeleton distance helper
		 */
		TsneRBFInterpolator(MocapReader * mr, SkeletonAnalizer * sa) : RBF2DInterpolator(mr, sa) {
		}

		virtual ~TsneRBFInterpolator() {}
		virtual Eigen::MatrixXd getMarkers(Eigen::VectorXi const & fids);
		virtual Eigen::MatrixXd getMerged2DProjection(MocapReader * sib);
};

#endif
