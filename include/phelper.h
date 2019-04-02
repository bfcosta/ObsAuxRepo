#ifndef PHELPER_H
#define PHELPER_H

#include <multireader.h>
#include <skanalizer.h>
#include <animation.h>

void writeBuffer(double * pvec, int size, Eigen::MatrixXd const & proj);

void readBuffer(double * pvec, int size, Eigen::MatrixXd & proj);

Eigen::VectorXi getControlPoints(Eigen::MatrixXd & proj);

void normalizeProjection(Eigen::MatrixXd & ret);

double phi_proj(double d2);

Eigen::MatrixXd getRbfProjection(MocapReader * mocap, SkeletonAnalizer * skdist, 
	Eigen::MatrixXd & pmat, Eigen::VectorXi & kfs);

Eigen::MatrixXd getRBFPMatrix(MocapReader * mocap, SkeletonAnalizer * skdist, 
	Eigen::MatrixXd & pmat, Eigen::VectorXi & kfs);
	
Eigen::MatrixXd getRBFP4UnitMocap(MultiReader * mocap, unsigned int idx, 
	SkeletonAnalizer * skdist, Eigen::MatrixXd & lmat, Eigen::MatrixXd & pmat, 
	Eigen::VectorXi & lkfs, Eigen::VectorXi & kfs);

Eigen::VectorXd getErrorVector(Eigen::MatrixXd const & proj, MocapReader * rd, SkeletonAnalizer * skdist);

Eigen::MatrixXd getGradientMatrix(Eigen::MatrixXd const & proj, MocapReader * rd, 
	SkeletonAnalizer * skdist, Eigen::VectorXd const & cerr, double alpha);

Eigen::MatrixXd improveProjection(Eigen::MatrixXd const & proj, MocapReader * rd, 
	SkeletonAnalizer * skdist, unsigned int iter = 200, double alpha = 0.5);

void loadAnimation(MocapReader * rd, Eigen::MatrixXd & ret);

extern Animation * skfanim;
#endif
