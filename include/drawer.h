#ifndef SHAPEDRAWER_H
#define SHAPEDRAWER_H

#include <GL/glu.h>

/**
 * \class ShapeDrawer
 *
 * \brief Useful solid drawing handler
 *
 * A base class for rendering solid shperes and cylinders with 
 * GLU interface. Useful for drawing skeleton, markers and controller.
 * 
 * @see http://www.glprogramming.com/red/chapter11.html
 * 
 */

class ShapeDrawer {

    protected:
		/// GLU quadric object
        GLUquadricObj * quadric;
        /// scope counter for memory management
        int * scope;

    public:
		/**
		 * creates a new object using a new GLU quadric 
		 * object for solid and smooth rendering.
		 */        
		ShapeDrawer() {
            quadric = gluNewQuadric();
            gluQuadricDrawStyle(quadric, GLU_FILL);
            gluQuadricNormals(quadric, GLU_SMOOTH);
            scope = new int(1);
		}
		
		/**
		 * \param that source object
		 * 
		 * creates a new object using a previous known GLU 
		 * quadric object for solid and smooth rendering.
		 */
        ShapeDrawer(const ShapeDrawer & that) {
            quadric = that.quadric;
            scope = that.scope;
             __sync_fetch_and_add(scope, 1);
            //(*scope)++;
        }
        
		/**
		 * \param that source object
		 * 
		 * assignment operator to use a previous known GLU 
		 * quadric object for solid and smooth rendering.
		 */
        virtual ShapeDrawer & operator=(const ShapeDrawer & that) {
            quadric = that.quadric;
            scope = that.scope;
             __sync_fetch_and_add(scope, 1);
            //(*scope)++;
            return *this;
        }

		/**
		 * Object destructor
		 */
        virtual ~ShapeDrawer () {
             __sync_fetch_and_sub(scope, 1);
            //(*scope)--;
            if (*scope ==0 ) {
                gluDeleteQuadric(quadric);
                delete scope;
            }
        }

		/**
		 * \param radius sphere radius
		 * \param slices sphere longitudinal divisions
		 * \param stacks sphere latitude divisions
		 * 
		 * draws a solid sphere centered at origin given current
		 * parameters
		 */
        void drawSphere(GLdouble radius, GLint slices, GLint stacks) {
            gluSphere(quadric , radius , slices , stacks);
        }
        
		/**
		 * \param radius cylinder radius in both extremes
		 * \param height cylinder length in Z axis
		 * \param slices angular divisions in XY plane
		 * \param stacks Z axis divisions
		 * 
		 * draws a solid cylinder along Z axis centered at origin 
		 * given current parameters.
		 */
        void drawCylinder(GLdouble radius, GLdouble height, GLint slices, GLint stacks) {
            gluCylinder(quadric, radius, radius, height, slices, stacks);
        }

		/**
		 * \param xbase length in X axis
		 * \param ybase length in Y axis
		 * \param height length in Z axis
		 * 
		 * draws a solid parallelepiped with XY rectangular shape along the Z axis and centered at origin.
		 */
        void drawParallelepiped(GLdouble xbase, GLdouble ybase, GLdouble height) {
			GLdouble x,y;
			x = xbase / 2;
			y = ybase / 2;
			
            glBegin(GL_QUADS);
            // XY plane lower
            glVertex3d( -x, -y, 0);
            glVertex3d( x, -y, 0);
            glVertex3d( x, y, 0);
            glVertex3d( -x, y, 0);
            // XY plane upper
            glVertex3d( -x, -y, height);
            glVertex3d( x, -y, height);
            glVertex3d( x, y, height);
            glVertex3d( -x, y, height);
            // XZ plane lower
            glVertex3d( -x, -y, 0);
            glVertex3d( x, -y, 0);
            glVertex3d( x, -y, height);
            glVertex3d( -x, -y, height);
            // XZ plane upper
            glVertex3d( -x, y, 0);
            glVertex3d( x, y, 0);
            glVertex3d( x, y, height);
            glVertex3d( -x, y, height);
            // YZ plane lower
            glVertex3d( -x, -y, 0);
            glVertex3d( -x, y, 0);
            glVertex3d( -x, y, height);
            glVertex3d( -x, -y, height);
            // YZ plane upper
            glVertex3d( x, -y, 0);
            glVertex3d( x, y, 0);
            glVertex3d( x, y, height);
            glVertex3d( x, -y, height);
            glEnd();
        }
};

#endif
