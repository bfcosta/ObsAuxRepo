#ifndef SKELETON_H
#define SKELETON_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <drawer.h>

/**
 * \class Skeleton
 *
 * \brief Describes the skeleton and its pose.
 *
 * A skeleton is a hierarchy of bones and joints. Here, joints are
 * modelled as translation matrices and bones are modeled as rotation
 * matrices. In an animation file, a joint is mostly static, although
 * some formats provide the ability to strech the bone length. Due to
 * this, most of the pose is described by the rotation matrices, known
 * here as bones.
 * 
 */

class Skeleton : public ShapeDrawer {

	protected:
	/// skeleton joint (translation) array
	std::vector<Eigen::Matrix4d> joint;
	/// skeleton bone (rotation) array
	std::vector<Eigen::Matrix4d> bone;
	/// a global rotation matrix for the skeleton root
	Eigen::Matrix4d orientation;
	/// name of the bone
	std::vector<std::string> name;
	/// index to parent joint
	std::vector<unsigned int> parent;
    /// skeleton color
    float scolor[4];
	/// bone thickness
	double thick;

	Eigen::Matrix4d getBoneMatrix(unsigned int i);

	public:
	Skeleton(double th = 1.0);
    //virtual ~Skeleton() {}
    void setColor(float r, float g, float b, float a = 1.0f);
    void setThickness(double th);
    //void addRoot(double tx, double ty, double tz, const std::string & n);
    //void addRoot(const Eigen::Matrix4d & offset, const std::string & n);
	void setRoot(double tx, double ty, double tz, const std::string & n);
	void setRoot(const Eigen::Matrix4d & offset, const std::string & n);
	void setRootPosition(const Eigen::Vector4d & r0);
	void setOrientation(double rx, double ry, double rz, const char * order);
	void setOrientation(const Eigen::Matrix4d & ori);
	void addJoint(double tx, double ty, double tz, 
		double rx, double ry, double rz, const char * order,
		unsigned int par, const std::string & n);
	void addJoint(const Eigen::Matrix4d & offset,
		const Eigen::Matrix4d & pRotation, unsigned int par, 
		const std::string & n);
	void draw(unsigned int glmode);
	Eigen::Vector4d getJointPosition(unsigned int i);
    Eigen::Matrix4d getRootOrientation();
    Eigen::Matrix4d getSkeletonRootTranslation();
	std::vector<unsigned int> getIndexHierarchy() {
		return parent;
	}
	std::vector<std::string> getBoneNames() {
		return name;
	}
	void setSkeleton(const Skeleton & s);
    void setPose(const Eigen::Matrix4d & jroot,
                 const Eigen::Matrix4d & ori, const std::vector<Eigen::Matrix4d> & brot);
	std::vector<unsigned int> getEndEffectorsIdx();
	std::vector<Eigen::Matrix4d> & getJointMatrices();
	std::vector<Eigen::Matrix4d> & getBoneMatrices();
	std::vector<Eigen::Matrix4d> getBonesTransfMatrices();
	std::vector<Eigen::Vector4d> getJointVectors();
	Eigen::Vector4d getCenter();
	double getRadius();
	void debug();
};

#endif
