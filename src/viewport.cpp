#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  Matrix3x3 translation = Matrix3x3::identity();
  translation(0, 2) = -((double)x - span);
  translation(1, 2) = -((double)y - span);
  Matrix3x3 scale = Matrix3x3::identity();
  scale(0, 0) = 1 / (span * 2);
  scale(1, 1) = 1 / (span * 2);
  set_canvas_to_norm(scale * translation);
  this->x = x;
  this->y = y;
  this->span = span; 

}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248
