/*
 * Implementations for Joint Based Characters.
 *
 * Started on October 29th, 2015 by Bryce Summers.
 */

#include "character.h"

#include "GL/glew.h"

namespace CMU462
{
   bool Joint :: calculateAngleGradient( Joint* goalJoint, Vector2D p, Vector2D q )
   {
      // TODO IMPLEMENT ME (TASK 2A)
       Vector2D u,diff,target;
       double gradient;
       bool influence=false;
       if(this==goalJoint)influence=true;
       for(int i=0;i<this->kids.size();i++)
       {
           if(this->kids.at(i)->calculateAngleGradient(goalJoint, p, q))
               influence=true;
       }
       if(influence)
       {
           target=p-currentCenter;//p是点击点（当前位置），currentCenter是父节点的中心坐标
           u=Vector2D(-target.y,target.x);//旋转90度
           diff=q-p;//q鼠标光标，
           gradient=dot(diff,u);
           ikAngleGradient=gradient;
           return true;
       }
       else{
           ikAngleGradient=0;
           return false;
       }

     // return false;
   }

   void Character :: reachForTarget( Joint* goalJoint,
                                     Vector2D sourcePoint,
                                     Vector2D targetPoint,
                                     double time )//指向点击点的指针，父节点源点位置（未转换），光标位置（目标点），当前动画时间
   {
      // TODO IMPLEMENT ME (TASK 2B)

       int step=10;
       for (int i=0;i<step;i++)
       {
           update(time);
           Vector3D v3=Vector3D(sourcePoint.x,sourcePoint.y,1);
           Vector3D mid=(goalJoint->currentTransformation)*v3;//将父节点源点转换成当前配置
           Vector2D v=Vector2D(mid.x,mid.y);
           root->calculateAngleGradient(goalJoint, v, targetPoint);
           for (int i=0;i<joints.size();i++)//遍历character的每一个子节点
           {
               double tau=0.000001;//时间步长
               double diff=tau*(joints.at(i)->ikAngleGradient);//移动的角度（时间乘以角度地变化率）
               double newAngle=joints.at(i)->getAngle(time)-diff;//新的角度
               joints.at(i)->setAngle(time, newAngle);//设置新的角度
           }
       }


   }


   void Joint :: integrate( double time, double timestep, Vector2D cumulativeAcceleration )
   {
      // TODO IMPLEMENT ME (TASK 3A)

       /***************task3A with preserving energy************/
//       double m,I;
//       const double a = 1;
       
//       Vector2D C;
//       physicalQuantities(m, I, C, center);
//       Vector2D mid=center-C;
//       double L=mid.norm();
//       double g=10;
//       double domega=-m*g*L*sin(theta)/I;//加速度
//       omega = (omega + domega*timestep)*exp(-a * timestep);//半隐式欧拉
//       theta+=omega*timestep;//半隐式欧拉
//       for(int i=0;i<kids.size();i++)
//       {
//           kids[i]->integrate(timestep, timestep, cumulativeAcceleration);
//       }
       /***************task3A with preserving energy************/

       double m, I;
       Vector2D C ;
       physicalQuantities(m, I, C , center);
       double L = (C  - center).norm()/100;
       double g = 3000;
       double psi = atan2(Vector2D(center - C ).x, -Vector2D(center - C ).y);//矢量Y与竖直方向的角度
       // using semi-implicit Euler
       // apply damping by multiplying exp(-a * timestep) where a = 1
       // also consider swing motion

       Vector2D y_perp=Vector2D(-Vector2D(center - C ).y,Vector2D(center - C ).x);//旋转中心到质心，逆时针旋转90度
       this->omega += -(dot(cumulativeAcceleration, y_perp) / (L*L) +
                        m * g * L * sin(this->theta-psi) / I) * timestep;
       this->theta += this->omega * timestep;
       this->omega *= exp(-timestep); //damping


       // Integrate the children.
       for( int i=0;i < kids.size();i++) {//求旋转引起的线性加速度

           // handle keyframed joints
           if( type == KEYFRAMED ) {

               // Lnew is different from L
               double sub_L=(center - (kids.at(i))->center).norm();
               double alpha = angle.evaluate(time);
               double alpha1 = angle.evaluate(time,1);
               double alpha2 = angle.evaluate(time,2);
               Vector3D accelJoint;
               accelJoint.x = L * (-alpha1*alpha1*sin(alpha)+alpha2*cos(alpha));//sin(alpha)的二阶导数
               accelJoint.y= L * (alpha1 * alpha1 * cos(alpha) + alpha2 * sin(alpha));
               accelJoint.z=0;
               accelJoint = currentTransformation * accelJoint;

               cumulativeAcceleration += Vector2D(accelJoint.x , accelJoint.y)*10;
           }

           (kids.at(i))->integrate(time, timestep, cumulativeAcceleration);
       }

       /***************task3B************/
   }

   void Character :: integrate( double time, double timestep )
   {
      // TODO IMPLEMENT ME (TASK 3B)
       Vector2D overallAccel=position.evaluate(time,2);//整体的加速度
       root->integrate(time, timestep, overallAccel);
      
   }

   void Character :: update( double time )
   {
      currentTransformation = Matrix3x3::translation( position.evaluate( time ) );
      root->update( time, currentTransformation );
   }

   void Character::draw( SVGRenderer* renderer, bool pick, Joint* hovered, Joint* selected )
   {
      root->draw( renderer, pick, hovered, selected );
   }

   // Every Joint is grouped with a circle representing the center.
   // The SVG file contains a hieracrhy of groups cooresponding to
   void Character::load_from_SVG(SVG & svg)
   {
      Group * root_group = static_cast<Group*>(svg.elements[0]);

      joints.clear();
      root = new Joint();
      root->index = joints.size();
      joints.push_back(root);

      root->parse_from_group(root_group, *this);
   }

   // The constructor sets the dynamic angle and velocity of
   // the joint to zero (at a perfect vertical with no motion)
   Joint :: Joint( void )
   : theta( 0. ), omega( 0. )
   {}

   Joint :: ~Joint( void )
   {
      // deallocate storage of SVG shapes
      for( int i = 0; i < shapes.size(); i++ )
      {
         delete shapes[i];
      }
   }

   double Joint::getAngle( double time )
   {
      if( type == DYNAMIC )
      {
         return theta;
      }

      // type == KEYFRAMED
      return angle( time );
   }

   void Joint::setAngle( double time, double value )
   {
      if( type == DYNAMIC )
      {
         theta = value;
         return;
      }

      // type == KEYFRAMED
      angle.setValue( time, value );
   }

   bool Joint::removeAngle(double time)
   {
      // Assuming times are on the integers only.
      return angle.removeKnot(time, .1);
   }

   void Joint::resetVelocity( void )
   {
      omega = 0.;
   }

   void Joint::resetDynamics( void )
   {
      theta = 0.;
      omega = 0.;
   }

   void Joint::update( double time, Matrix3x3 parentTransformation )
   {
      // Calculate the cumulative transformation by composing the
      // transformation of the parent with a rotation around the
      // joint center.

      double alpha = getAngle( time );
      if( type == DYNAMIC )
      {
         // A pendulum should hang straight down; its angle should
         // not be affected by the rotation of any joints above it.
         alpha -= parentTransformation.getRotation();
      }

      Matrix3x3 R = Matrix3x3::rotation( alpha );
      Matrix3x3 A = Matrix3x3::translation(  center );
      Matrix3x3 B = Matrix3x3::translation( -center );
      currentParentTransformation = parentTransformation;
      currentTransformation = parentTransformation * A * R * B;

      Vector3D c( center.x, center.y, 1. );
      c = currentTransformation * c;
      currentCenter = Vector2D( c.x/c.z, c.y/c.z );

      for( vector<Joint*>::iterator joint  = kids.begin(); joint != kids.end(); joint ++ )
      {
         (*joint)->update( time, currentTransformation );
      }
   }

   inline double mod1(double in) { return in > 1.0 ? in - 1 : in;}


   void Joint::draw( SVGRenderer* renderer, bool pick, Joint* hovered, Joint* selected )
   {
      for( vector<SVGElement*>::iterator shape = shapes.begin(); shape != shapes.end(); shape++ )
      {
         // make a copy of the original style for this shape,
         // so that we can restore it if it's changed for either
         // picking or hovering/selection
         Style originalStyle = (*shape)->style;

         if( pick )
         {
            // compute and apply the picking color
            Color pickColor = Color::fromPickIndex( index+1 );
            (*shape)->style.strokeColor = pickColor;
            (*shape)->style.fillColor   = pickColor;
         }
         else if( this == selected )
         {
            // highlight the fill and stroke color
            // Selected joints are drawn with a contrasting color.
            Color& fillColor( (*shape)->style.fillColor );
            fillColor.r = mod1(fillColor.r + .5);
            fillColor.g = mod1(fillColor.g + .5);
            fillColor.b = mod1(fillColor.b + .5);


            Color& strokeColor( (*shape)->style.fillColor );
            /*
               fillColor.r = 1.;
               fillColor.g = 1.;
               fillColor.b = 1.;
               */
         }
         else if( this == hovered )
         {
            // highlight the fill color
            Color& fillColor( (*shape)->style.fillColor );
            fillColor.r += .35;
            fillColor.g += .35;
            fillColor.b += .35;
         }

         renderer->pushTransformation();
         renderer->concatenateTransformation( currentTransformation );
         renderer->draw_element( *shape );
         renderer->popTransformation();

         // restore the original style
         (*shape)->style = originalStyle;
      }

      for( vector<Joint*>::iterator joint = kids.begin(); joint != kids.end(); joint ++ )
      {
         (*joint)->draw( renderer, pick, hovered, selected );
      }
   }

   void Joint::parse_from_group(Group * G, Character & C)
   {
      // Attempt to parse the starting group containing the current joint's data.
      vector<SVGElement*> & elements = G -> elements;
      vector<SVGElement*>::iterator iter = elements.begin();

      // Parse this particular group's data.
      Group * joint_group = dynamic_cast<Group*>(*iter);

      // Base Case, this is childless joint if this is not a group.
      if(joint_group == NULL)
      {
         // one shape for each element,
         //excluding the final circle which just specifies the joint
         int nShapes = elements.size()-1;
         shapes.resize( nShapes );
         for( int i = 0; i < nShapes; i++ )
         {
            shapes[i] = (*iter)->copy();
            iter++;
         }
         Circle * center_circle =  static_cast<Circle*>(*iter);
         center = center_circle -> center;
         setJointType( center_circle );
         return;
      }

      // one shape for each element,
      // excluding the final circle which just specifies the joint.
      int nShapes = joint_group->elements.size()-1; 
      shapes.resize( nShapes );
      for( int i = 0; i < nShapes; i++ )
      {
         shapes[i] = ((joint_group->elements)[i])->copy();
      }
      SVGElement * svg_circle = (joint_group->elements)[nShapes];
      Circle * center_circle  = dynamic_cast<Circle*>(svg_circle);

      if(center_circle == NULL)
      {
         cout << svg_circle << endl;
         cout << svg_circle->type << endl;
         cerr << "ERROR: The Center of rotation circle was not found during joint parsing.\n";
         exit(0);
      }

      center = center_circle -> center;
      setJointType( center_circle );

      iter++;


      // Parse the sub joints.
      while(iter != elements.end())
      {
         Joint * child = new Joint();

         // Allocate new joints for the children in the Character's joints vector.
         child->index = C.joints.size();
         C.joints.push_back(child);

         // Link each of these joints to this joint's kids array.
         kids.push_back(child);

         Group * child_group = static_cast<Group*>(*iter);

         // Parse the child.
         child -> parse_from_group(child_group, C);

         // Continue on to the next child.
         iter++;
      }

      // Done parsing subjoints.

   }

   void Joint :: setJointType( Circle* circle )
   {
      Color& c( circle->style.fillColor );

      type = KEYFRAMED;

      if( c.r == 1. && c.g == 1. && c.b == 1. )
      {
         type = DYNAMIC;
      }
   }

   void Joint :: physicalQuantities( double& m, // mass
                                     double& I, // moment of inertia
                                     Vector2D& c, // centroid
                                     Vector2D center
                                     ) const
   {
      m = 0.;
      I = 0.;
      c = Vector2D( 0., 0. );

      for( int k = 0; k < shapes.size(); k++ )
      {
         double   mk = shapes[k]->mass();
         double   Ik = shapes[k]->momentOfInertia( center );
         Vector2D ck = shapes[k]->centroid();

         m += mk;
         I += Ik;
         c += mk*ck;
      }

      c /= m;
   }
}
