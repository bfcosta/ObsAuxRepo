#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <kfselector.h>
#include <bvhreader.h>
#include <phelper.h>
#include <emscripten/emscripten.h>

MocapReader * mocap = NULL;
Animation * skfanim = NULL;
MultiReader * multi = NULL;
std::vector<MocapReader *> mpool;
Eigen::MatrixXd mproj, rproj;
int bframe = 0, nframes = 0;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief load BVH contents in database.
 * @param argv string with BVH file contents.
 * @return number of frames in BVH
 */
int add2DB(char * argv) {
    // save bvh file
    std::string input(argv);
    std::ofstream out("saved.bvh", std::ofstream::trunc);
    out << input;
    out.close();
	MocapReader * rd = new BVHReader();
	rd->readFile("saved.bvh");
	mpool.push_back(rd);
	return rd->getNumberOfFrames();	
}

/**
 * @brief clear database contents.
 * @return zero if database exists; one otherwise
 */
int clearDB() {
	mpool.clear();
	if (multi != NULL) {
		delete multi;
		multi = NULL;
		return 0;
	}
	return 1;
}

/**
 * @brief select keyframes in database using ROLS+MDS.
 * @param pvout output buffer to write frames (id+projection)
 * @param nkfs number of keyframes
 * @param nsamples use maximum number of samples in DB to select keyframes
 * @return number of keyframes
 */
int computeKFs(double * pvout, unsigned int nkfs, unsigned int nsamples) {
	
	// build multireader
	if (nkfs > nsamples) nsamples = nkfs;
	if (multi == NULL) {
		multi = new MultiReader(nsamples);
		//if (mocap == NULL) mocap = multi;
	}
	for (unsigned int i=multi->getNumberOfReaders(); i < mpool.size(); ++i) {
		multi->addMocap(mpool.at(i));
	}
	multi->reassign();
	//find keyframes in multireader
	std::cout << "requested # of KFs = " << nkfs << std::endl;
	double pc = multi->getNumberOfFrames();
	std::cout << "Sample size = " << pc << std::endl;
	std::cout << "Maximum sample size = " << nsamples << std::endl;
	pc = nkfs / pc;
	std::cout << "KF % = " << pc*100 << std::endl;
	RolsRbfpEucMdsRbfSelector rolsrbfpmds(multi, pc);
	rolsrbfpmds.printStats();
	// flush results
	Eigen::VectorXi kfs = rolsrbfpmds.getKeyframeVector();
	Eigen::MatrixXd proj2 = rolsrbfpmds.getFramesProjection();
	mproj = proj2.block(0,0,proj2.rows(),3);
	proj2 = Eigen::ArrayXXd::Zero(kfs.size(), mproj.cols());
	for (unsigned int i=0; i < kfs.size(); ++i) {
		mproj(kfs(i),2) = 1;
		proj2(i,0) = multi->getFrameIdBySampleId(kfs(i));
		proj2(i,1) = mproj(kfs(i), 0);
		proj2(i,2) = mproj(kfs(i), 1);
	}
	//std::cout << "keyframes # = " << kfs.size() << std::endl;
	//std::cout << "buffer address = " << pvout << std::endl;
	//std::cout << "keyframes projection" << std::endl;
	//std::cout << proj2 << std::endl;
	writeBuffer(pvout, 3*nkfs, proj2);
	loadAnimation(multi, mproj);
	return kfs.size();
}

/**
 * @brief read keyframes and related projections to build a new SKF interpolation
 * @param pvin projection buffer
 * @param nkfs number of elements in buffer
 * @return number of keyframes read
 */
int setKFs(double * pvin, unsigned int nkfs) {
	unsigned int fid;
	
	Eigen::MatrixXd proj = Eigen::ArrayXXd::Zero(nkfs, 3);
	Eigen::VectorXi kfs = Eigen::ArrayXi::Zero(nkfs);
	readBuffer(pvin, 3*nkfs, proj);
	//~ std::cout << "last projection" << std::endl;
	//~ for (unsigned int i=0; i < mproj.rows(); ++i)
		//~ if (mproj(i,2) != 0)
			//~ std::cout << i << ": " << mproj.row(i) << std::endl;
	//~ std::cout << "keyframes projection" << std::endl;
	//~ std::cout << proj << std::endl;
	// get keyframes in sample
	for (unsigned int i=0; i < kfs.size(); ++i) {
		fid = proj(i,0);
		kfs(i) = multi->getSampleIdByFrameId(fid);
	}
	//~ std::cout << "keyframes = " << kfs.transpose() << std::endl;
	// clear keyframes
	for (unsigned int i=0; i < mproj.rows(); ++i) 
		mproj(i,2) = 0;
	// new keyframes
	for (unsigned int i=0; i < kfs.size(); ++i) {
		mproj(kfs(i),0) = proj(i,1);
		mproj(kfs(i),1) = proj(i,2);
		mproj(kfs(i),2) = 1;
	}
	// load animation
	//~ std::cout << "current projection" << std::endl;
	//~ for (unsigned int i=0; i < mproj.rows(); ++i)
		//~ if (mproj(i,2) != 0)
			//~ std::cout << i << ": " << mproj.row(i) << std::endl;
	loadAnimation(multi, mproj);
	return kfs.size();
}

/**
 * @brief get RBF Projection
 * @param pvout buffer with RBFP: a three-column row-wise matrix
 * where third column is the control point flag
 * @param base first frame of the projection
 * @param nf total of frames to be projected
 * @return number of projected frames
 */
int getRBFP(double * pvout, int base, int nf) {
	GeodesicSkeletonAnalizer skan;
	Eigen::MatrixXd tmp;
	Eigen::VectorXi lctps;
	MocapReader * rd;
	unsigned int pf = 0, ff = 0, lf = 0;
	
	bframe = base;
	nframes = nf;
	Eigen::VectorXi ctps = getControlPoints(mproj);
	rproj = Eigen::ArrayXXd::Zero(nf,mproj.cols());
	Eigen::MatrixXd lmat = getRBFPMatrix(multi, &skan, mproj, ctps);
	for (unsigned int i=0; i < multi->getNumberOfReaders(); ++i) {
		ff = lf;
		rd = multi->getReader(i);
		lf = ff + rd->getNumberOfFrames();
		if (((base >= ff) && (base < lf)) || ((base < ff) && (pf < nf))) {
			tmp = getRBFP4UnitMocap(multi, i, &skan, lmat, mproj, lctps, ctps);
			for (unsigned int j=0; j < tmp.rows(); ++j) 
				if (((ff + j) >= base) && (pf < nf))
					rproj.row(pf++) = tmp.row(j);
		}
	}
	writeBuffer(pvout, 3*nf, rproj);
	return pf;
}

/**
 * @brief get optimized RBF Projection
 * @param pvout buffer with RBFP: a three-column row-wise matrix
 * where third column is the control point flag. Size of buffer was given
 * in last call of getRBFP.
 * @return number of projected frames
 */
int getRBFPOpt(double * pvout) {
	GeodesicSkeletonAnalizer skan;
	Eigen::MatrixXd tmp0, tmp1;
	MocapReader * rd;
	unsigned int pf1 = 0, pf0 = 0, ff = 0, lf = 0;
	
	int base = bframe;
	int nf = nframes;
	for (unsigned int i=0; i < multi->getNumberOfReaders(); ++i) {
		ff = lf;
		rd = multi->getReader(i);
		lf = ff + rd->getNumberOfFrames();
		if (((base >= ff) && (base < lf)) || ((base < ff) && (pf0 < nf))) {
			tmp0 = Eigen::ArrayXXd::Zero(lf - ff, rproj.cols()); 
			for (unsigned int j=0; j < tmp0.rows(); ++j) 
				if (((ff + j) >= base) && (pf0 < nf))
					tmp0.row(j) = rproj.row(pf0++); // rproj was set on last getRBFP call
			tmp1 = improveProjection(tmp0, rd, &skan);
		}
	}
	writeBuffer(pvout, 3*nf, tmp1);
	return pf1;
}

/**
 * @brief get error measure between two skeleton poses
 * @param qvin1 quaternion buffer for first pose
 * @param qvin2 quaternion buffer for second pose
 * @param qsize quaternion buffer size
 * @return pose distance
 */
double getPoseError(double * qvin1, double * qvin2, int qsize) {
	GeodesicSkeletonAnalizer skan;
	Skeleton s0, s1;
	Eigen::Quaterniond q;
	int cnt = 0;

	// rotations
	if (multi == NULL) {
		std::cout << "No mocap was loaded." << std::endl;
		return -1;
	}
	//std::cout << "Building skeletons." << std::endl;
	s0 = multi->getSkeleton();
	s1 = multi->getSkeleton();
	Eigen::Matrix4d or0 = s0.getRootOrientation();
	Eigen::Matrix4d root0 = s0.getJointMatrices().at(0);
	std::vector<std::string> bnames = s0.getBoneNames();
	std::vector<Eigen::Matrix4d> bmat0 = s0.getBoneMatrices();
	std::vector<Eigen::Matrix4d> bmat1 = s1.getBoneMatrices();
	//std::cout << "Copying skeleton poses." << std::endl;
	for (unsigned int i=0; (i < bnames.size()) && (cnt < qsize); ++i) {
		if (bnames.at(i).compare("leaf") != 0) {
			q.x() = qvin1[cnt];
			q.y() = qvin1[cnt+1];
			q.z() = qvin1[cnt+2];
			q.w() = qvin1[cnt+3];
			bmat0.at(i).block<3,3>(0,0) = q.normalized().toRotationMatrix();
			q.x() = qvin2[cnt];
			q.y() = qvin2[cnt+1];
			q.z() = qvin2[cnt+2];
			q.w() = qvin2[cnt+3];
			bmat1.at(i).block<3,3>(0,0) = q.normalized().toRotationMatrix();
			cnt += 4;
		}
	}
	if (cnt != qsize) {
		std::cout << "Quaternions and rotations matrices mismatch." << std::endl;
		std::cout << "Number of quaternions = " << qsize/4 << std::endl;
		std::cout << "Quaternions read = " << cnt/4 << std::endl;
		std::cout << "Number of rotation matrices = " << bmat0.size() << std::endl;
	}
	s0.setPose(root0, or0, bmat0);
	s1.setPose(root0, or0, bmat1);
	return (skan.getPoseDistance(s0, s1) / multi->getNumberOfDOFs());
}

/**
 * @brief get skeleton pose for a given 2D coordinate for the last loaded spatial keyframe animation
 * @param qvout output skeleton rotation buffer
 * @param qsize size of output skeleton rotation buffer
 * @param xpos X coordinate
 * @param ypos Y coordinate
 * @return number of rotations degrees of freedom
 */
int unprojectPose(double * qvout, int qsize, double xpos, double ypos) {
	Eigen::Vector4d v;
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
	return qvec.size();	
}

#ifdef __cplusplus
}
#endif
