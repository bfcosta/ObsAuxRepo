#include <iostream>
#include <fstream>
#include <string>
#include <phelper.h>

/**
 * @brief write contents of a given matrix into a raw buffer
 * @param pvec buffer to be written
 * @param size buffer size
 * @param proj source matrix
 */
void writeBuffer(double * pvec, int size, Eigen::MatrixXd const & proj) {
	int tmp;
	for (unsigned int i=0; i < proj.rows(); ++i)
		for (unsigned int j=0; j < proj.cols(); ++j) {
			tmp = i*proj.cols() + j;
			if (tmp < size)
				pvec[tmp] = proj(i,j);
			else
				break;
		}
}

/**
 * @brief read contents of a raw buffer and write them into a matrix
 * @param pvec buffer to be read
 * @param size buffer size
 * @param proj target matrix
 */
void readBuffer(double * pvec, int size, Eigen::MatrixXd & proj) {
	int tmp;
	for (unsigned int i=0; i < proj.rows(); ++i) 
		for (unsigned int j=0; j < proj.cols(); ++j) {
			tmp = i*proj.cols() + j;
			if (tmp < size) 
				proj(i,j) = pvec[tmp];
			else
				break;
		}
}

/**
 * @brief get control points
 * @param proj projected data rowwise matrix
 * @return keyframe vector
 */
Eigen::VectorXi getControlPoints(Eigen::MatrixXd & proj) {

    Eigen::VectorXi sample = Eigen::ArrayXi::Zero(proj.col(2).sum());
    unsigned int idx = 0;
    for (unsigned int i=0; i < proj.rows(); ++i)
		if (proj(i,2) != 0)
			sample(idx++) = i;
	return sample;
}

/**
 * @brief normalize projection to square [0, 1]
 * @param ret projection matrix to be normalized
 */
void normalizeProjection(Eigen::MatrixXd & ret) {
	double min0, max0;
	double min1, max1;
	// normalize projection coordinates
	min0 = max0 = ret(0,0);
	min1 = max1 = ret(0,1);
	for (unsigned int i=1; i < ret.rows(); ++i) {
		if (min0 > ret(i,0)) min0 = ret(i,0);
		else if (max0 < ret(i,0)) max0 = ret(i,0);
		if (min1 > ret(i,1)) min1 = ret(i,1);
		else if (max1 < ret(i,1)) max1 = ret(i,1);
	}
	for (unsigned int i=0; i < ret.rows(); ++i) {
		ret(i,0) = (ret(i,0) - min0) / (max0 - min0);
		ret(i,1) = (ret(i,1) - min1) / (max1 - min1);
	}
}

/**
 * @brief return the implemented radial basis function kernel as quadratic multiquadrics
 * with c and epsilon equals one.
 * @param d2 distance as input
 * @return function output
 */
double phi_proj(double d2) {
	//return std::sqrt(1 + d2*d2);  // multiquadrics - elisa
	return (1 + std::pow(d2, 2.5));
}

/**
 * @brief get RBF projection onto a 2D plane given the control points 
 * and their corresponding projections
 * @param mocap mocap session handler
 * @param skdist pose distance definition a SkeletonAnalizer object
 * @param pmat projection matrix
 * @param kfs control points vector
 * @return a projection matrix with all poses
 */
Eigen::MatrixXd getRbfProjection(MocapReader * mocap, SkeletonAnalizer * skdist, 
	Eigen::MatrixXd & pmat, Eigen::VectorXi & kfs) {
	
	Eigen::MatrixXd proj = pmat;
	int frames = mocap->getNumberOfFrames();
	// find lambda
	Eigen::MatrixXd fmat = Eigen::ArrayXXd::Zero(kfs.size(), kfs.size());
	Eigen::MatrixXd ymat = Eigen::ArrayXXd::Zero(kfs.size(), 2);
	Skeleton sk0, sk1;
	double dist;
	
	for (int i=0; i < fmat.rows(); ++i)
		for (int j=i+1; j < fmat.cols(); ++j) {
			sk0 = mocap->getFrameSkeleton(kfs(i));
			sk1 = mocap->getFrameSkeleton(kfs(j));
			dist = std::sqrt(skdist->getPoseDistance(sk0, sk1));
			fmat(j,i) = fmat(i,j) = phi_proj(dist);
		}
	for (int i=0; i < ymat.rows(); ++i) {
		ymat(i,0) = proj(kfs(i),0);
		ymat(i,1) = proj(kfs(i),1);
	}
	Eigen::MatrixXd lmat = fmat.inverse() * ymat;
	
	// project remainder frames by RBF
	Eigen::VectorXd rbfvec = Eigen::ArrayXd::Zero(kfs.size());
	int c = 0;
	for (int i=0; i < frames; ++i) {
		if (i == kfs(c)) {
			if (c != (kfs.size() - 1)) ++c;
		} else {
		//if (proj(i,2) == 0) {
			for (int j=0; j < rbfvec.size(); ++j) {
				sk0 = mocap->getFrameSkeleton(i);
				sk1 = mocap->getFrameSkeleton(kfs(j));
				dist = std::sqrt(skdist->getPoseDistance(sk0, sk1));
				rbfvec(j) = phi_proj(dist);
			}
			proj.row(i).head(2) = rbfvec.transpose() * lmat;
		}
	}
	normalizeProjection(proj);
	//std::cout << "RBF projection" << std::endl << proj << std::endl;
	return proj;
}
	
/**
 * @brief get RBFP matrix that outlines the corresponding projections
 * @param mocap mocap session handler
 * @param skdist pose distance definition a SkeletonAnalizer object
 * @param pmat seed projection matrix
 * @param kfs control points vector
 * @return RBFP matrix
 */
Eigen::MatrixXd getRBFPMatrix(MocapReader * mocap, SkeletonAnalizer * skdist, 
	Eigen::MatrixXd & pmat, Eigen::VectorXi & kfs) {
	
	Eigen::MatrixXd proj = pmat;
	int frames = mocap->getNumberOfFrames();
	// find lambda
	Eigen::MatrixXd fmat = Eigen::ArrayXXd::Zero(kfs.size(), kfs.size());
	Eigen::MatrixXd ymat = Eigen::ArrayXXd::Zero(kfs.size(), 2);
	Skeleton sk0, sk1;
	double dist;
	
	for (int i=0; i < fmat.rows(); ++i)
		for (int j=i+1; j < fmat.cols(); ++j) {
			sk0 = mocap->getFrameSkeleton(kfs(i));
			sk1 = mocap->getFrameSkeleton(kfs(j));
			dist = std::sqrt(skdist->getPoseDistance(sk0, sk1));
			fmat(j,i) = fmat(i,j) = phi_proj(dist);
		}
	for (int i=0; i < ymat.rows(); ++i) {
		ymat(i,0) = proj(kfs(i),0);
		ymat(i,1) = proj(kfs(i),1);
	}
	Eigen::MatrixXd lmat = fmat.inverse() * ymat;
	return lmat;	
}

/**
 * @brief get RBF projection for a unit mocap handler
 * @param mocap unit mocap session handler
 * @param idx index of unit mocap handler
 * @param skdist pose distance definition a SkeletonAnalizer object
 * @param lmat RBFP matrix
 * @param pmat seed projection matrix
 * @param lkfs control points for chosen unit mocap handler
 * @param kfs control points vector
 * @return projection matrix
 */
Eigen::MatrixXd getRBFP4UnitMocap(MultiReader * mocap, unsigned int idx,
	SkeletonAnalizer * skdist, Eigen::MatrixXd & lmat, 
	Eigen::MatrixXd & pmat, Eigen::VectorXi & lkfs, Eigen::VectorXi & kfs) {
	
	Skeleton sk0, sk1;
	double dist;
	MocapReader * rd = mocap->getReader(idx);
	int frames = rd->getNumberOfFrames();
	Eigen::VectorXi fsample = mocap->getFrameSampleByUnitReader(idx);
	unsigned int fbase = mocap->getFrameBaseByUnitReader(idx);
	Eigen::MatrixXd proj = Eigen::ArrayXXd::Zero(rd->getNumberOfFrames(), pmat.cols());
	
	// find global keyframes in local unit mocap range
	unsigned int lf = 0;
	for (unsigned int i=0; i < kfs.size(); ++i) {
		if ((kfs(i) >= fbase) && (kfs(i) < (fbase + fsample.size()))) {
			++lf;
		}
	}
	lkfs = Eigen::ArrayXi::Zero(lf);
	lf = 0;
	for (unsigned int i=0; i < kfs.size(); ++i) {
		if ((kfs(i) >= fbase) && (kfs(i) < (fbase + fsample.size()))) {
			lkfs(lf) = fsample(kfs(i) - fbase);
			proj.row(lkfs(lf)) = pmat.row(kfs(i));
			++lf;			
		}
	}

	// project remainder frames by RBF
	Eigen::VectorXd rbfvec = Eigen::ArrayXd::Zero(kfs.size());
	int c = 0;
	for (int i=0; i < frames; ++i) {
		if (i == lkfs(c)) {
			if (c != (lkfs.size() - 1)) ++c;
		} else {
			for (int j=0; j < rbfvec.size(); ++j) {
				sk0 = rd->getFrameSkeleton(i);
				sk1 = mocap->getFrameSkeleton(kfs(j));
				dist = std::sqrt(skdist->getPoseDistance(sk0, sk1));
				rbfvec(j) = phi_proj(dist);
			}
			proj.row(i).head(2) = rbfvec.transpose() * lmat;
		}
	}
	//normalizeProjection(proj);
	return proj;
}

/**
 * @brief get the error vector of a reconstructed RBF/SKF animation for a given
 * frame projection matrix
 * @param proj frame projection matrix
 * @param mocap multi session handler
 * @param rd mocap session handler
 * @param skdist pose distance definition a SkeletonAnalizer object
 * @return error vector
 */
Eigen::VectorXd getErrorVector(Eigen::MatrixXd const & proj, MocapReader * rd, SkeletonAnalizer * skdist) {

	Eigen::VectorXd cerr = Eigen::ArrayXd::Zero(proj.rows());
	Eigen::VectorXd v0;
	Eigen::Vector4d cpos;
	Skeleton sk0, sk1;
	
	cpos(2) = 0;
	for (unsigned int i=0; i < cerr.size(); ++i) {
		sk0 = rd->getFrameSkeleton(i);
		v0 = sk0.getJointPosition(0);
		// current error
		cpos.head(2) = proj.row(i).head(2);
		cpos(3) = 1;
		skfanim->setController(cpos);
		skfanim->animate();
		sk1 = skfanim->getSkeleton();
		sk1.setRootPosition(v0);
		cerr(i) = skdist->getPoseDistance(sk0, sk1);
	}
	return cerr;
}

/**
 * @brief get the gradient matrix of a 2D pose projection with inverse projection
 * given by RBF/SKF.
 * @param proj pose projection matrix
 * @param rd mocap session handler
 * @param skdist pose distance definition a SkeletonAnalizer object
 * @param cerr error vector
 * @param alpha gain parameter
 * @return gradient matrix
 */
Eigen::MatrixXd getGradientMatrix(Eigen::MatrixXd const & proj, MocapReader * rd, 
	SkeletonAnalizer * skdist, Eigen::VectorXd const & cerr, double alpha) {
		
	Eigen::MatrixXd error = Eigen::ArrayXXd::Zero(proj.rows(), 5);
	Eigen::MatrixXd grad = Eigen::ArrayXXd::Zero(proj.rows(), 2);
	Skeleton sk0, sk1;
	Eigen::VectorXd v0;
	Eigen::Vector4d cpos, cbpos, capos;

	// get NSEW delta error (just for first)
	cpos(2) = 0;
	capos(2) = 0;
	cbpos(2) = 0;
	sk0 = rd->getFrameSkeleton(0);
	v0 = sk0.getJointPosition(0);
	cpos.head(2) = proj.row(0).head(2);
	cpos(3) = 1;
	capos.head(2) = proj.row(1).head(2);
	capos(3) = 1;
	capos -= cpos;
	error(0,4) = capos.norm(); // circle radius
	capos = cpos;
	capos(0) += alpha*error(0,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(0,0) = skdist->getPoseDistance(sk0, sk1);
	error(0,0) -= cerr(0);
	// delta west error
	capos = cpos;
	capos(0) -= alpha*error(0,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(0,1) = skdist->getPoseDistance(sk0, sk1);
	error(0,1) -= cerr(0);
	// delta north error
	capos = cpos;
	capos(1) += alpha*error(0,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(0,2) = skdist->getPoseDistance(sk0, sk1);
	error(0,2) -= cerr(0);
	// delta south error
	capos = cpos;
	capos(1) -= alpha*error(0,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(0,3) = skdist->getPoseDistance(sk0, sk1);
	error(0,3) -= cerr(0);
	for (unsigned int i=1; i < error.rows()-1; ++i) {
		// delta east error
		sk0 = rd->getFrameSkeleton(i);
		v0 = sk0.getJointPosition(0);
		cpos.head(2) = proj.row(i).head(2);
		cpos(3) = 1;
		cbpos.head(2) = proj.row(i-1).head(2);
		cbpos(3) = 1;
		capos.head(2) = proj.row(i+1).head(2);
		capos(3) = 1;
		cbpos -= cpos;
		capos -= cpos;
		error(i,4) = (cbpos.norm() + capos.norm()) / 2; // circle radius
		capos = cpos;
		capos(0) += alpha*error(i,4);
		skfanim->setController(capos);
		skfanim->animate();
		sk1 = skfanim->getSkeleton();
		sk1.setRootPosition(v0);
		error(i,0) = skdist->getPoseDistance(sk0, sk1);
		error(i,0) -= cerr(i);
		// delta west error
		capos = cpos;
		capos(0) -= alpha*error(i,4);
		skfanim->setController(capos);
		skfanim->animate();
		sk1 = skfanim->getSkeleton();
		sk1.setRootPosition(v0);
		error(i,1) = skdist->getPoseDistance(sk0, sk1);
		error(i,1) -= cerr(i);
		// delta north error
		capos = cpos;
		capos(1) += alpha*error(i,4);
		skfanim->setController(capos);
		skfanim->animate();
		sk1 = skfanim->getSkeleton();
		sk1.setRootPosition(v0);
		error(i,2) = skdist->getPoseDistance(sk0, sk1);
		error(i,2) -= cerr(i);
		// delta south error
		capos = cpos;
		capos(1) -= alpha*error(i,4);
		skfanim->setController(capos);
		skfanim->animate();
		sk1 = skfanim->getSkeleton();
		sk1.setRootPosition(v0);
		error(i,3) = skdist->getPoseDistance(sk0, sk1);
		error(i,3) -= cerr(i);
	}
	// last
	sk0 = rd->getFrameSkeleton(error.rows()-1);
	v0 = sk0.getJointPosition(0);
	cpos.head(2) = proj.row(error.rows()-1).head(2);
	cpos(3) = 1;
	capos.head(2) = proj.row(error.rows()-2).head(2);
	capos(3) = 1;
	capos -= cpos;
	error(error.rows()-1,4) = capos.norm(); // circle radius
	capos = cpos;
	capos(0) += alpha*error(error.rows()-1,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(error.rows()-1,0) = skdist->getPoseDistance(sk0, sk1);
	error(error.rows()-1,0) -= cerr(error.rows()-1);
	// delta west error
	capos = cpos;
	capos(0) -= alpha*error(error.rows()-1,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(error.rows()-1,1) = skdist->getPoseDistance(sk0, sk1);
	error(error.rows()-1,1) -= cerr(error.rows()-1);
	// delta north error
	capos = cpos;
	capos(1) += alpha*error(error.rows()-1,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(error.rows()-1,2) = skdist->getPoseDistance(sk0, sk1);
	error(error.rows()-1,2) -= cerr(error.rows()-1);
	// delta south error
	capos = cpos;
	capos(1) -= alpha*error(error.rows()-1,4);
	skfanim->setController(capos);
	skfanim->animate();
	sk1 = skfanim->getSkeleton();
	sk1.setRootPosition(v0);
	error(error.rows()-1,3) = skdist->getPoseDistance(sk0, sk1);
	error(error.rows()-1,3) -= cerr(error.rows()-1);
	for (unsigned int i=0; i < error.rows(); ++i) {
		grad(i,0) = 0;
		grad(i,1) = 0;
		if ((error(i,0) < 0) || (error(i,1) < 0)) {
			grad(i,0) = error(i,1) - error(i,0);
		}
		if ((error(i,2) < 0) || (error(i,3) < 0)) {
			grad(i,1) = error(i,3) - error(i,2);
		}
	}
	for (unsigned int i=0; i < grad.rows(); ++i) {
		if (proj(i,2) == 1)
			grad.row(i) *= 0;
		if ((grad(i,0) != 0) || (grad(i,1) != 0)) {
			grad.row(i).normalize();
			grad.row(i) *= error(i,4);
		}
	}
	return grad;
}

/**
 * @brief improve a given projection matrix
 * @param proj pose projection matrix
 * @param rd mocap session handler
 * @param skdist pose distance definition a SkeletonAnalizer object
 * @param iter number of gradient descent iterations
 * @param alpha gradient descent gain parameter
 * @return improved frame projection matrix
 */
Eigen::MatrixXd improveProjection(Eigen::MatrixXd const & proj, MocapReader * rd, 
	SkeletonAnalizer * skdist, unsigned int iter, double alpha) {
		
	Eigen::MatrixXd grad;
	Eigen::VectorXd cerr;
	Eigen::MatrixXd ret = proj;
	double tini, total0, total1;
	unsigned int i;
	
	cerr = getErrorVector(ret, rd, skdist);
	tini = total0 = cerr.sum();
	//~ std::cout << "original projection error = " << total0 << std::endl;
	for (i=0; i < iter; ++i) {
		grad = getGradientMatrix(ret, rd, skdist, cerr, alpha);
		ret.block(0, 0, ret.rows(), grad.cols()) += alpha * grad;
		cerr = getErrorVector(ret, rd, skdist);	
		total1 = cerr.sum();
		//~ std::cout << "delta projection error = " << total1 << std::endl;
		if (total1 == total0) break;
		total0 = total1;
	}
	total1 = 100 * (tini - total0) / tini;
	tini /= rd->getNumberOfDOFs() * rd->getNumberOfFrames();
	total0 /= rd->getNumberOfDOFs() * rd->getNumberOfFrames();
	//~ std::cout << "original projection error = " << tini << std::endl;
	//~ std::cout << "delta projection error = " << total0 << std::endl;
	std::cout << "total of iterations = " << i << std::endl;
	//~ std::cout << "error gain = " << total1 << " %" << std::endl;
	//normalizeProjection(ret);
	return ret;
}

/**
 * @brief load animation given a projection matrix
 * @param rd mocap session handler
 * @param ret projection matrix
 */
void loadAnimation(MocapReader * rd, Eigen::MatrixXd & ret) {
		
	Eigen::Vector4d v;
	Skeleton s;

	// load animation
	s = rd->getFrameSkeleton(0);
	if (skfanim != NULL)
		delete skfanim;
	skfanim = new Animation();
	skfanim->setSkeleton(s);
	v << 0,0,0,1;
	for (unsigned int i=0; i < ret.rows(); ++i) {
		if (ret(i,2) != 0) {
			v(0) = ret(i,0);
			v(1) = ret(i,1);
			s = rd->getFrameSkeleton(i);	
			skfanim->addNewKeyFrame(v, s);
		}
	}
	v(0) = ret(0,0);
	v(1) = ret(0,1);
	skfanim->setController(v);
	skfanim->setContext();
}
