#include <Eigen/Dense>
#include <vector>

/**
 * \class RbfBaseSolver
 *
 * \brief common data and methods for rbf solvers
 *
 * Base class for solving RBF interpolations. Provides useful data
 * and methods for its child classes. The RBF implementation
 * interpolates an input set of dimension 2 or 3 into a value. 
 * Default input set dimension is 2. The input set elements are
 * also known as markers, representing points in plane or space
 * 
 * \note This class shouldn't be used directly, as no output is provided.
 * Use RbfSimpleInterp or RbfSkelInterp as you like.
 *
 */

class RbfBaseSolver {
	protected:
	/// data input set for RBF
    std::vector<Eigen::VectorXd> markers;
	/// matrix which provides the interpolation
    Eigen::MatrixXd inverse;
	/// input set dimension
	unsigned int din;
	/// flag to avoid redundant RBF computation
    bool redundant;
    /// compact support radius
    double csradius;
    //QDebug dbgout;

    double phi(Eigen::VectorXd const & v1, Eigen::VectorXd const & v2);
    virtual void addRedundant();
    virtual void delRedundant();
    double getValue(Eigen::VectorXd const & sol, Eigen::VectorXd const & p);

	public:
    RbfBaseSolver(unsigned int d = 2) {
        din = d; 
        redundant = false;
        csradius = 1;
    }

	/**
	 *  clears all data
	 */
	virtual ~RbfBaseSolver() {
		reset();
	}

    void addMarker(Eigen::VectorXd const & m);
    void addMarker(std::vector<Eigen::VectorXd> const & mvec);
    void addMarker(std::vector<Eigen::Vector4d> const & mvec);
	void prepare();

	/**
	 *  clears all data
	 */
	virtual void reset() {
        din = 2; 
        redundant = false;
		markers.clear();
	}

	/**
	 *  reset data input dimension to 3
	 */
	void set3DimEntrySet() {
		din = 3;
	}

	/**
	 *  set compact support radius
	 * \param r new compact support radius
	 */
    void setCsradius(double r) {
		csradius = r;
	}
    Eigen::VectorXd solve(Eigen::VectorXd const & m);
    Eigen::VectorXd getDistVec(Eigen::VectorXd const & p);    
	virtual void dump();
};

/**
 * \class RbfSimpleInterp
 *
 * \brief Models and solves RBF interpolations of 1D output
 *
 * A class for solving RBF interpolations of 1D output. 
 * If provided one image value to each given marker, it calculates
 * RBF coefficients which will be used to interpolate a new value.
 * If output data is constant, no interpolation is calculated
 * and the constant value is returned.
 * 
 * \note Useful mostly for debugging purposes. Not used on
 * Spatial KF.
 */ 
class RbfSimpleInterp: public RbfBaseSolver {
	protected:
	///vector of marker outputs, one for each given marker
    std::vector<double> image;
	/// calculated rbf polynomial coefficients
    Eigen::VectorXd rbfCoefs;
	/// a flag to avoid calculation, if output data is constant
    bool constval;
    
    virtual void addRedundant();
    virtual void delRedundant();

	public:
	RbfSimpleInterp(unsigned int d = 2): RbfBaseSolver(d) {}

	/**
	 *  clears all data
	 */
	virtual ~RbfSimpleInterp() {
		reset();
	}

    void addImage(double im);
    void addImage(std::vector<double> const & imvec);
	void solve();
	void resolve(Eigen::VectorXd const & im);
    double getImage(Eigen::VectorXd const & m);

	/**
	 *  clears all data
	 */
	virtual void reset() {
		this->RbfBaseSolver::reset();
		image.clear();
	}
	virtual void dump();
};

/**
 * \class RbfSkelInterp
 *
 * \brief Models and solves RBF interpolations of N dim output
 *
 * A class for solving RBF interpolations of N dim output. 
 * If provided N image values to each given marker, it calculates
 * RBF coefficients which will be used to interpolate N new value.
 * If output data is constant, no interpolation is calculated
 * and the constant value is returned. On Spatial KF, N=12 mostly.
 * 
 * \note this class implements Spatial KF core calculus.
 */
class RbfSkelInterp: public RbfBaseSolver {
	protected:
	/// vector of outputs, one vector of rotation matrix for each given marker
    std::vector< std::vector<Eigen::Matrix4d> > skel;
	/// calculated rbf polynomial coefficients for each rotation matrix
    std::vector< std::vector<Eigen::VectorXd> > rbfBoneCoefs;
	/// a flag to avoid calculation, if output data is constant
    std::vector< std::vector<bool> > constel;

    virtual void addRedundant();
    virtual void delRedundant();
    void orthonormalize(Eigen::Matrix4d & n);

	public:

	RbfSkelInterp(): RbfBaseSolver() {}

	/**
	 *  clears all data
	 */
	virtual ~RbfSkelInterp() {
		reset();
	}
    void addSkel(std::vector<Eigen::Matrix4d> const & rvec);
	void solve();
    std::vector<Eigen::Matrix4d> const getNewSkel(Eigen::Vector4d const & m);

	/**
	 *  clears all data
	 */
	virtual void reset() {
		this->RbfBaseSolver::reset();
		skel.clear();
		rbfBoneCoefs.clear();
		constel.clear();
	}
	virtual void dump();
};
