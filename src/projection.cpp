#include <iostream>
#include <fstream>
#include <string>
#include <kfselector.h>
#include <bvhreader.h>
#include <animation.h>
#include <skanalizer.h>
#include <emscripten/emscripten.h>

MocapReader * mocap = NULL;
Animation * skfanim = NULL;

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

    //~ Eigen::VectorXi sample = Eigen::ArrayXi::Zero(proj.rows());
    //~ for (unsigned int i=0; i < proj.rows(); ++i)
		//~ sample(i) = proj(i,0);
	//~ return sample;
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
 * @brief get the error vector of a reconstructed RBF/SKF animation for a given
 * frame projection matrix
 * @param proj frame projection matrix
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
	SkeletonAnalizer * skdist, unsigned int iter = 200, double alpha = 0.5) {
		
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
	normalizeProjection(ret);
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


#ifdef __cplusplus
extern "C" {
#endif


//~ //int EMSCRIPTEN_KEEPALIVE saveBvh(char * argv) {
//~ int saveBvh(char * argv) {
    //~ std::string input(argv);
    //~ std::ofstream out("saved.bvh", std::ofstream::trunc);
    //~ out << input;
    //~ out.close();
	//~ return 0;
//~ }

/**
 * @brief load BVH contents, save handler for future use and gets total of frames.
 * @param argv string with BVH file contents.
 * @return number of frames in BVH
 */
unsigned int getNumberOfFrames(char * argv) {
    // save bvh file
    std::string input(argv);
    std::ofstream out("saved.bvh", std::ofstream::trunc);
    out << input;
    out.close();
    // get projection
	if (mocap != NULL) {
		delete mocap;
		mocap = NULL;
	}
	mocap = new BVHReader();
	mocap->readFile("saved.bvh");
	return mocap->getNumberOfFrames();
}

/**
 * @brief get keyframes and MDS projection for last loaded BVH
 * @param pc percentage of frames. Expected range is (0; 100)
 * @param pvec output buffer vector in which projection is loaded as a row wise matrix
 * @param size size of buffer vector elements
 * @return address of buffer vector
 */
double * getBvhProjection(double pc, double * pvec, int size) {
	if (mocap == NULL) {
		for (int i=0; i < size; ++i)
			pvec[i] = -1;
	} else {
		RolsRbfpEucMdsRbfSelector rolsrbfpmds(mocap, pc/100);
		rolsrbfpmds.printStats();
		Eigen::MatrixXd proj2 = rolsrbfpmds.getFramesProjection();
		Eigen::MatrixXd proj = proj2.block(0,0,proj2.rows(),3);
		Eigen::VectorXi kfs = rolsrbfpmds.getKeyframeVector();
		//Eigen::VectorXi kfs = rolsrbfpmds.getKeyframes(unsigned int nkf);
		for (unsigned int i=0; i < kfs.size(); ++i)
			proj(kfs(i),2) = 1;
		writeBuffer(pvec, size, proj);
		loadAnimation(mocap, proj);
	}
	return pvec;
}

/**
 * @brief get RBFP given a seed projection with control points
 * @param pvin buffer with seed projection: a three-column row-wise matrix
 * where third column is the control point flag
 * @param pvout buffer with RBFP: a three-column row-wise matrix
 * where third column is the control point flag
 * @param size buffer size
 * @return output buffer address
 */
double * getRbfSmoothedProjection(double * pvin, double * pvout, int size) {
	GeodesicSkeletonAnalizer skan;
	
	Eigen::MatrixXd proj = Eigen::ArrayXXd::Zero(size/3, 3);
	readBuffer(pvin, size, proj);
	loadAnimation(mocap, proj);
	Eigen::VectorXi ctps = getControlPoints(proj);
	Eigen::MatrixXd rproj = getRbfProjection(mocap, &skan, proj, ctps);
	writeBuffer(pvout, size, rproj);
	return pvout;
}

/**
 * @brief get forward and backward improved projection with control points.
 * Mocap handler has to have been built with getNumberOfFrames. Also, control points 
 * must be chosen by getBvhProjection.
 * @param pvin buffer with out-of-box projection
 * @param pvout buffer with improved projection
 * @param size buffer size
 * @return output buffer address
 */
double * getImprovedProjection(double * pvin, double * pvout, int size) {
	GeodesicSkeletonAnalizer skan;
	
	Eigen::MatrixXd proj = Eigen::ArrayXXd::Zero(size/3, 3);
	readBuffer(pvin, size, proj);
	loadAnimation(mocap, proj);
	Eigen::MatrixXd iproj = improveProjection(proj, mocap, &skan);
	writeBuffer(pvout, size, iproj);
	return pvout;
}

/**
 * @brief get skeleton pose error between a given frame and a pose built by
 * SKF interpolation given a controller position in 2D.
 * @param pvin buffer with a given projection: a three-column row-wise matrix
 * where third column is the control point flag
 * @param size pvin buffer size
 * @param xpos X position on the controller plane
 * @param ypos Y position on the controller plane
 * @param frame mocap frame id
 * @return calculated error between frame id and SKF interpolation
 */
double getPoseError(double * pvin, int size, double xpos, double ypos, int frame) {
	GeodesicSkeletonAnalizer skan;
	Eigen::Vector4d v;
	Eigen::VectorXd v0;
	Skeleton s0, s1;
	Animation anim;
	
	Eigen::MatrixXd proj = Eigen::ArrayXXd::Zero(size/3, 3);
	readBuffer(pvin, size, proj);
	// load animation
	s1 = mocap->getFrameSkeleton(frame);
	v0 = s1.getJointPosition(0);
	anim.setSkeleton(s1);
	v << 0,0,0,1;
	for (unsigned int i=0; i < proj.rows(); ++i) {
		if (proj(i,2) != 0) {
			v(0) = proj(i,0);
			v(1) = proj(i,1);
			s0 = mocap->getFrameSkeleton(i);	
			anim.addNewKeyFrame(v, s0);
		}
	}
	v(0) = xpos;
	v(1) = ypos;
	anim.setController(v);
	anim.setContext();
	anim.animate();
	s0 = anim.getSkeleton();
	s0.setRootPosition(v0);
	return (skan.getPoseDistance(s0, s1) / mocap->getNumberOfDOFs());
}

/**
 * @brief get skeleton pose error between a given frame and a pose built by
 * the last loaded SKF interpolation given a controller position in 2D.
 * @param xpos X position on the controller plane
 * @param ypos Y position on the controller plane
 * @param frame mocap frame id
 * @return calculated error between frame id and SKF interpolation
 */
double getLoadedPoseError(double xpos, double ypos, int frame) {
	GeodesicSkeletonAnalizer skan;
	Eigen::Vector4d v;
	Eigen::VectorXd v0;
	Skeleton s0, s1;
	
	// load animation
	s1 = mocap->getFrameSkeleton(frame);
	v << 0,0,0,1;
	v(0) = xpos;
	v(1) = ypos;
	if (skfanim != NULL) {
		skfanim->setController(v);
		skfanim->setContext();
		skfanim->animate();
		s0 = skfanim->getSkeleton();
	} else {
		return -1;
	}
	v0 = s1.getJointPosition(0);
	s0.setRootPosition(v0);
	return (skan.getPoseDistance(s0, s1) / mocap->getNumberOfDOFs());
}

/**
 * @brief get skeleton pose for a given 2D coordinate
 * @param pvin input projection buffer
 * @param isize size of input projection buffer
 * @param qvout output skeleton rotation buffer
 * @param qsize size of output skeleton rotation buffer
 * @param jvout output skeleton position buffer
 * @param jsize size of output skeleton position buffer
 * @param xpos X coordinate
 * @param ypos Y coordinate
 * @return number of rotations degrees of freedom
 */
int getSkfSkeleton(double * pvin, int isize, double * qvout, int qsize, 
	double * jvout, int jsize, double xpos, double ypos) {
		
	Eigen::Vector4d v;
	Eigen::VectorXd v0;
	Skeleton s0;
	Animation anim;
	Eigen::Quaterniond quat;
	int bones = 0, bcnt = 0;

	Eigen::MatrixXd proj = Eigen::ArrayXXd::Zero(isize/3, 3);
	readBuffer(pvin, isize, proj);
	// load animation
	s0 = mocap->getFrameSkeleton(0);
	v0 = s0.getJointPosition(0);
	anim.setSkeleton(s0);
	v << 0,0,0,1;
	for (unsigned int i=0; i < proj.rows(); ++i) {
		if (proj(i,2) != 0) {
			v(0) = proj(i,0);
			v(1) = proj(i,1);
			s0 = mocap->getFrameSkeleton(i);	
			anim.addNewKeyFrame(v, s0);
		}
	}
	v(0) = xpos;
	v(1) = ypos;
	anim.setController(v);
	anim.setContext();
	anim.animate();
	s0 = anim.getSkeleton();
	// rotations
	std::vector<Eigen::Matrix4d> mvec = s0.getBoneMatrices();
	std::vector<std::string> bnames = s0.getBoneNames();
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) 
			++bones;
	}
	Eigen::VectorXd qvec = Eigen::ArrayXd::Zero(bones * 4);
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			quat = mvec.at(i).block<3,3>(0,0);
			qvec(4*bcnt) = quat.x();
			qvec(4*bcnt + 1) = quat.y();
			qvec(4*bcnt + 2) = quat.z();
			qvec(4*bcnt + 3) = quat.w();
			++bcnt;
		}
	}
	if (qsize > qvec.size()) 
		qsize = qvec.size();
	for (int i=0; i < qsize; ++i)
		qvout[i] = qvec(i);
	// joints
	std::vector<Eigen::Vector4d> jvec = s0.getJointVectors();
	Eigen::VectorXd posvec = Eigen::ArrayXd::Zero(bones * 3);
	bcnt = 0;
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			posvec(3*bcnt) = jvec.at(i)(0);
			posvec(3*bcnt + 1) = jvec.at(i)(1);
			posvec(3*bcnt + 2) = jvec.at(i)(2);
			++bcnt;
		}
	}
	if (jsize > posvec.size()) 
		jsize = posvec.size();
	for (int i=0; i < jsize; ++i)
		jvout[i] = posvec(i);
	return qvec.size();
}

/**
 * @brief get skeleton pose for a given 2D coordinate for the last loaded spatial keyframe animation
 * @param qvout output skeleton rotation buffer
 * @param qsize size of output skeleton rotation buffer
 * @param jvout output skeleton position buffer
 * @param jsize size of output skeleton position buffer
 * @param xpos X coordinate
 * @param ypos Y coordinate
 * @return number of rotations degrees of freedom
 */
int getLoadedSkfSkeleton(double * qvout, int qsize, double * jvout, int jsize, double xpos, double ypos) {
		
	Eigen::Vector4d v;
	Eigen::VectorXd v0;
	Skeleton s0;
	Eigen::Quaterniond quat;
	int bones = 0, bcnt = 0;

	v << 0,0,0,1;
	v(0) = xpos;
	v(1) = ypos;
	if (skfanim != NULL) {
		skfanim->setController(v);
		skfanim->setContext();
		skfanim->animate();
		s0 = skfanim->getSkeleton();
	} else {
		return 0;
	}
	// rotations
	std::vector<Eigen::Matrix4d> mvec = s0.getBoneMatrices();
	std::vector<std::string> bnames = s0.getBoneNames();
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) 
			++bones;
	}
	Eigen::VectorXd qvec = Eigen::ArrayXd::Zero(bones * 4);
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			quat = mvec.at(i).block<3,3>(0,0);
			qvec(4*bcnt) = quat.x();
			qvec(4*bcnt + 1) = quat.y();
			qvec(4*bcnt + 2) = quat.z();
			qvec(4*bcnt + 3) = quat.w();
			++bcnt;
		}
	}
	if (qsize > qvec.size()) 
		qsize = qvec.size();
	for (int i=0; i < qsize; ++i)
		qvout[i] = qvec(i);
	// joints
	std::vector<Eigen::Vector4d> jvec = s0.getJointVectors();
	Eigen::VectorXd posvec = Eigen::ArrayXd::Zero(bones * 3);
	bcnt = 0;
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			posvec(3*bcnt) = jvec.at(i)(0);
			posvec(3*bcnt + 1) = jvec.at(i)(1);
			posvec(3*bcnt + 2) = jvec.at(i)(2);
			++bcnt;
		}
	}
	if (jsize > posvec.size()) 
		jsize = posvec.size();
	for (int i=0; i < jsize; ++i)
		jvout[i] = posvec(i);
	return qvec.size();
}

/**
 * @brief get mocap skeleton given by a frame id
 * @param qvout output buffer of skeleton rotations
 * @param qsize size of output buffer of skeleton rotations
 * @param jvout output buffer of skeleton positions
 * @param jsize size of output buffer of skeleton positions
 * @param frameid frame id
 * @return number of rotations degrees of freedom
 */
int getMocapSkeleton(double * qvout, int qsize, double * jvout, int jsize, int frameid) {
	Skeleton s0;
	Eigen::Quaterniond quat;
	int bones = 0, bcnt = 0;

	s0 = mocap->getFrameSkeleton(frameid);
	std::vector<std::string> bnames = s0.getBoneNames();
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) 
			++bones;
	}
	// rotations
	std::vector<Eigen::Matrix4d> mvec = s0.getBoneMatrices();
	Eigen::VectorXd qvec = Eigen::ArrayXd::Zero(bones * 4);
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			quat = mvec.at(i).block<3,3>(0,0);
			qvec(4*bcnt) = quat.x();
			qvec(4*bcnt + 1) = quat.y();
			qvec(4*bcnt + 2) = quat.z();
			qvec(4*bcnt + 3) = quat.w();
			++bcnt;
		}
	}
	if (qsize > qvec.size()) 
		qsize = qvec.size();
	for (int i=0; i < qsize; ++i)
		qvout[i] = qvec(i);
	// joints
	std::vector<Eigen::Vector4d> jvec = s0.getJointVectors();
	Eigen::VectorXd posvec = Eigen::ArrayXd::Zero(bones * 3);
	bcnt = 0;
	for (unsigned int i=0; i < bnames.size(); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			posvec(3*bcnt) = jvec.at(i)(0);
			posvec(3*bcnt + 1) = jvec.at(i)(1);
			posvec(3*bcnt + 2) = jvec.at(i)(2);
			++bcnt;
		}
	}
	if (jsize > posvec.size()) 
		jsize = posvec.size();
	for (int i=0; i < jsize; ++i)
		jvout[i] = posvec(i);
	return qvec.size();
}

#ifdef __cplusplus
}
#endif
