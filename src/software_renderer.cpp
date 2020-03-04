#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <time.h> 
#include <stdio.h> 
#include <chrono>
#include <ctime>
#include "triangulation.h"
#include <cassert>
#include <array>
using namespace std;

namespace CS248 {

  static std::chrono::duration<double> elapsed_seconds1;
  static std::chrono::duration<double> elapsed_seconds2;
  static std::chrono::duration<double> elapsed_secondsC;
  static std::chrono::duration<double> elapsed_secondsCPP;
  static std::chrono::duration<double> elapsed_seconds5;
  // Implements SoftwareRenderer //

  // fill a sample location with color
  void SoftwareRendererImp::fill_sample(int sx, int sy, const Color& color) {
    size_t start = 4 * (sx + sy * target_w);
    render_target[start] = (uint8_t)(color.r * 255);
    render_target[start + 1] = (uint8_t)(color.g * 255);
    render_target[start + 2] = (uint8_t)(color.b * 255);
    render_target[start + 3] = (uint8_t)(color.a * 255);
  }

  //sx = sx * 4; sy = sy * target_w * 4;
  void SoftwareRendererImp::fill_sample2(int sx, int sy, const Color& color) {
    int start = sx + sy;
    render_target[start] = (uint8_t)(color.r * 255);
    render_target[start + 1] = (uint8_t)(color.g * 255);
    render_target[start + 2] = (uint8_t)(color.b * 255);
    render_target[start + 3] = (uint8_t)(color.a * 255);
  }

  // fill samples in the entire pixel specified by pixel coordinates
  void SoftwareRendererImp::fill_pixel(int x, int y, const Color& color) {

    // Task 2: Re-implement this function

    // check bounds
    if (x < 0 || x >= target_w) return;
    if (y < 0 || y >= target_h) return;

    Color pixel_color;
    float inv255 = 1.0f / 255.0f;
    pixel_color.r = render_target[4 * (x + y * target_w)] * inv255;
    pixel_color.g = render_target[4 * (x + y * target_w) + 1] * inv255;
    pixel_color.b = render_target[4 * (x + y * target_w) + 2] * inv255;
    pixel_color.a = render_target[4 * (x + y * target_w) + 3] * inv255;

    pixel_color = ref->alpha_blending_helper(pixel_color, color);

    render_target[4 * (x + y * target_w)] = (uint8_t)(pixel_color.r * 255);
    render_target[4 * (x + y * target_w) + 1] = (uint8_t)(pixel_color.g * 255);
    render_target[4 * (x + y * target_w) + 2] = (uint8_t)(pixel_color.b * 255);
    render_target[4 * (x + y * target_w) + 3] = (uint8_t)(pixel_color.a * 255);

  }

  void SoftwareRendererImp::draw_svg(SVG& svg) {

    // set top level transformation
    transformation = canvas_to_screen;

    // draw all elements
    for (size_t i = 0; i < svg.elements.size(); ++i) {
      draw_element(svg.elements[i]);
    }

    // draw canvas outline
    Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
    Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
    Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
    Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

    rasterize_line((float)a.x, (float)a.y, (float)b.x, (float)b.y, Color::Black);
    rasterize_line((float)a.x, (float)a.y, (float)c.x, (float)c.y, Color::Black);
    rasterize_line((float)d.x, (float)d.y, (float)b.x, (float)b.y, Color::Black);
    rasterize_line((float)d.x, (float)d.y, (float)c.x, (float)c.y, Color::Black);

    // resolve and send to render target
    resolve();

  }

  void SoftwareRendererImp::set_sample_rate(size_t sample_rate) {

    // Task 2: 
    // You may want to modify this for supersampling support
    this->sample_rate = sample_rate;

  }

  void SoftwareRendererImp::set_render_target(unsigned char* render_target,
    size_t width, size_t height) {

    // Task 2: 
    // You may want to modify this for supersampling support
    this->render_target = render_target;
    this->target_w = width;
    this->target_h = height;

  }

  void SoftwareRendererImp::draw_element(SVGElement* element) {

    // Task 3 (part 1):
    // Modify this to implement the transformation stack

    switch (element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
    }

  }


  // Primitive Drawing //

  void SoftwareRendererImp::draw_point(Point& point) {

    Vector2D p = transform(point.position);
    rasterize_point((float)p.x, (float)p.y, point.style.fillColor);

  }

  void SoftwareRendererImp::draw_line(Line& line) {

    Vector2D p0 = transform(line.from);
    Vector2D p1 = transform(line.to);
    rasterize_line((float)p0.x, (float)p0.y, (float)p1.x, (float)p1.y, line.style.strokeColor);

  }

  void SoftwareRendererImp::draw_polyline(Polyline& polyline) {

    Color c = polyline.style.strokeColor;

    if (c.a != 0) {
      int nPoints = (int)polyline.points.size();
      for (int i = 0; i < nPoints - 1; i++) {
        Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
        Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
        rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
      }
    }
  }

  void SoftwareRendererImp::draw_rect(Rect& rect) {

    Color c;

    // draw as two triangles
    float x = rect.position.x;
    float y = rect.position.y;
    float w = rect.dimension.x;
    float h = rect.dimension.y;

    Vector2D p0 = transform(Vector2D(x, y));
    Vector2D p1 = transform(Vector2D(x + w, y));
    Vector2D p2 = transform(Vector2D(x, y + h));
    Vector2D p3 = transform(Vector2D(x + w, y + h));

    // draw fill
    c = rect.style.fillColor;
    if (c.a != 0) {
      rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
      rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
    }

    // draw outline
    c = rect.style.strokeColor;
    if (c.a != 0) {
      rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
      rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
      rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
      rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
    }

  }

  void SoftwareRendererImp::draw_polygon(Polygon& polygon) {

    Color c;

    // draw fill
    c = polygon.style.fillColor;
    if (c.a != 0) {

      // triangulate
      vector<Vector2D> triangles;
      triangulate(polygon, triangles);

      // draw as triangles
      for (size_t i = 0; i < triangles.size(); i += 3) {
        Vector2D p0 = transform(triangles[i + 0]);
        Vector2D p1 = transform(triangles[i + 1]);
        Vector2D p2 = transform(triangles[i + 2]);
        rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
      }
    }

    // draw outline
    c = polygon.style.strokeColor;
    if (c.a != 0) {
      int nPoints = polygon.points.size();
      for (int i = 0; i < nPoints; i++) {
        Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
        Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
        rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
      }
    }
  }

  void SoftwareRendererImp::draw_ellipse(Ellipse& ellipse) {

    // Extra credit 

  }

  void SoftwareRendererImp::draw_image(Image& image) {

    Vector2D p0 = transform(image.position);
    Vector2D p1 = transform(image.position + image.dimension);

    rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
  }

  void SoftwareRendererImp::draw_group(Group& group) {

    for (size_t i = 0; i < group.elements.size(); ++i) {
      draw_element(group.elements[i]);
    }

  }

  // Rasterization //

  // The input arguments in the rasterization functions 
  // below are all defined in screen space coordinates

  void SoftwareRendererImp::rasterize_point(float x, float y, Color color) {

    // fill in the nearest pixel
    //int sx = (int)floor(x);
    //int sy = (int)floor(y);
    int sx = (int)x;
    int sy = (int)y;
    // check bounds
    if (sx < 0 || sx >= target_w) return;
    if (sy < 0 || sy >= target_h) return;

    // fill sample - NOT doing alpha blending!
    // TODO: Call fill_pixel here to run alpha blending
    fill_sample(sx, sy, color);
  }

  void SoftwareRendererImp::rasterize_line1(float x0f, float y0f,
    float x1f, float y1f,
    Color color) {
    int x0 = (int)x0f;
    int y0 = (int)y0f;
    int x1 = (int)x1f;
    int y1 = (int)y1f;
    int dy = y1 - y0;
    int dx = x1 - x0;
    float t = 0.0f;
    if (abs(dx) > abs(dy)) {
      float m = (float)dy / (float)dx;
      t = y0;
      dx = (dx < 0) ? -1 : 1;
      m *= dx;
      while (x0 != x1) {
        fill_sample(x0, (int)t, color);
        x0 += dx;
        t += m;
      }
    }
    else {
      float m = (float)dx / (float)dy;
      t = x0;
      dy = (dy < 0) ? -1 : 1;
      m *= dy;
      while (y0 != y1) {
        fill_sample((int)t, y0, color);
        y0 += dy;
        t += m;
      }
    }
  }

  void SoftwareRendererImp::rasterize_line2(float x0f, float y0f,
    float x1f, float y1f,
    Color color) {
    int x0 = (int)x0f;
    int y0 = (int)y0f;
    int x1 = (int)x1f;
    int y1 = (int)y1f;
    int dy = y1 - y0;
    int dx = x1 - x0;
    int stepx = 4;
    int stepy = target_w << 2;
    if (dy < 0) {
      dy = -dy;
      stepy = -stepy;
    }
    if (dx < 0) {
      dx = -dx;
      stepx = -stepx;
    }

    x0 <<= 2;
    x1 <<= 2;
    y0 *= target_w * 4;
    y1 *= target_w * 4;

    //int start = 4 * (sx + sy * target_w);
    if (dx > dy) {
      int fraction = -dx;
      while (x0 != x1) {
        fill_sample2(x0, y0, color);
        //int advance = fraction / dx;

        if (fraction >= 0) {
          y0 += stepy;
          fraction -= dx;
        }
        x0 += stepx;
        fraction += dy;
      }
    }
    else {
      int fraction = -dy;
      while (y0 != y1) {
        fill_sample2(x0, y0, color);
        //int advance = fraction / dy;

        if (fraction >= 0) {
          x0 += stepx;
          fraction -= dy;
        }
        y0 += stepy;
        fraction += dx;
      }
    }
  }

  void SoftwareRendererImp::rasterize_lineC(int x0, int y0, int x1, int y1, Color color) {

    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = (dx > dy ? dx : -dy) / 2, e2;

    for (;;) {
      fill_sample(x0, y0, color);
      if (x0 == x1 && y0 == y1) break;
      e2 = err;
      if (e2 > -dx) { err -= dy; x0 += sx; }
      if (e2 < dy) { err += dx; y0 += sy; }
    }
  }

  void SoftwareRendererImp::rasterize_lineCPP(int x1, int y1, int x2, int y2, Color color) {

    // Bresenham's line algorithm
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if (steep)
    {
      std::swap(x1, y1);
      std::swap(x2, y2);
    }

    if (x1 > x2)
    {
      std::swap(x1, x2);
      std::swap(y1, y2);
    }

    const float dx = x2 - x1;
    const float dy = fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for (int x = (int)x1; x <= maxX; x++)
    {
      if (steep)
      {
        fill_sample(y, x, color);
      }
      else
      {
        fill_sample(x, y, color);
      }

      error -= dy;
      if (error < 0)
      {
        y += ystep;
        error += dx;
      }
    }
  }

  typedef
    union
  {
    int32_t i;
    struct
    {
      int16_t lo; // endian-specific!
      int16_t hi;
    };
  } fixed_point;

  void SoftwareRendererImp::rasterize_line5(int x1, int y1, int x2, int y2, Color color) {

    int dy = y2 - y1;
    int dx = x2 - x1;

    fixed_point t = { 0 };
    if (abs(dx) > abs(dy)) {
      int32_t m = ((int32_t)dy << 16) / dx; // slope as fixed point
      t.i = y1 << 16;
      dx = (dx < 0) ? -1 : 1;
      m *= dx;
      while (x1 != x2) {
        fill_sample(x1, t.hi, color);
        x1 += dx;
        t.i += m;
      }
    }
    else {
      if (dy == 0) return;
      int32_t m = ((int32_t)dx << 16) / dy; // slope as fixed point
      t.i = x1 << 16;
      dy = (dy < 0) ? -1 : 1;
      m *= dy;
      while (y1 != y2) {
        fill_sample(t.hi, y1, color);
        y1 += dy;
        t.i += m;
      }
    }
  }

  void SoftwareRendererImp::rasterize_line(float x0f, float y0f,
    float x1f, float y1f,
    Color color) {

    //     short count = 1000000;
    // 
    //     auto start = ::std::chrono::system_clock::now();
    // 
    //     for (int i = 0; i < count; i++) {
    //         rasterize_line1(x0f, y0f, x1f, y1f, color);
    //     }
    //     auto end = std::chrono::system_clock::now();
    //     elapsed_seconds1 += (end - start);
    // 
    //     start = std::chrono::system_clock::now();
    //     for (int i = 0; i < count; i++) {
    //         rasterize_line2(x0f, y0f, x1f, y1f, color);
    //     }
    //     end = std::chrono::system_clock::now();
    //     elapsed_seconds2 += (end - start);
    //     
    //     start = std::chrono::system_clock::now();
    //     for (int i = 0; i < count; i++) {
    //         rasterize_lineC(x0f, y0f, x1f, y1f, color);
    //     }
    //     end = std::chrono::system_clock::now();
    //     elapsed_secondsC += (end - start);
    // 
    //     start = std::chrono::system_clock::now();
    //     for (int i = 0; i < count; i++) {
    //         rasterize_lineCPP(x0f, y0f, x1f, y1f, color);
    //     }
    //     end = std::chrono::system_clock::now();
    //     elapsed_secondsCPP += (end - start);
    // 
    //     start = std::chrono::system_clock::now();
    //     for (int i = 0; i < count; i++) {
    rasterize_line5(x0f, y0f, x1f, y1f, color);
    //     }
    //     end = std::chrono::system_clock::now();
    //     elapsed_seconds5 += (end - start);
    // 
    //     ::std::cout << "elapsed time1: " << elapsed_seconds1.count() << endl;
    //     ::std::cout << "elapsed time2: " << elapsed_seconds2.count() << endl;
    //     ::std::cout << "elapsed timeC: " << elapsed_secondsC.count() << endl;
    //     ::std::cout << "elapsed timeCPP: " << elapsed_secondsCPP.count() << endl;
    //     ::std::cout << "elapsed time5: " << elapsed_seconds5.count() << endl;
  }

  struct Point2D {
    float x;
    float y;
  };

  static bool cmp(const Point2D& a, const Point2D& b)
  {
    return (a.y < b.y || (a.y == b.y && a.x > b.x));
  }

  static float round_down(float v)//round down to the nearest sample point
  {
    return static_cast<int>(v - 0.5) + 0.5;
  }

  void SoftwareRendererImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {

    Vector2D v0((double)x1 - x0, (double)y1 - y0);
    Vector2D v1((double)x2 - x1, (double)y2 - y1);

    if (cross(v0, v1) < 0) {
      v0 = -v0;
      swap(x0, x1);
      swap(y0, y1);
    }
    Vector2D v2 = -(v1 + v0);
    auto xbound = std::minmax({ x0,x1,x2 });
    auto ybound = std::minmax({ y0,y1,y2 });
    float starty = static_cast<int>(ybound.first) + 0.5;
    float startx = static_cast<int>(xbound.first) + 0.5;
    float A0 = v0.y;
    float B0 = -v0.x;
    float C0 = y0 * v0.x - x0 * v0.y;

    float A1 = v1.y;
    float B1 = -v1.x;
    float C1 = y1 * v1.x - x1 * v1.y;

    float A2 = v2.y;
    float B2 = -v2.x;
    float C2 = y2 * v2.x - x2 * v2.y;

    float L0 = A0 * startx + B0 * starty + C0;
    float L1 = A1 * startx + B1 * starty + C1;
    float L2 = A2 * startx + B2 * starty + C2;
    for (float y = starty; y <= round_down(ybound.second); y += 1) {
      float L0_y = L0;
      float L1_y = L1;
      float L2_y = L2;
      for (float x = startx; x <= round_down(xbound.second); x += 1) {
        if (L0_y <= 0 && L1_y <= 0 && L2_y <= 0) {
          fill_sample(x, y, color);
        }
        L0_y += A0;
        L1_y += A1;
        L2_y += A2;
      }
      L0 += B0;
      L1 += B1;
      L2 += B2;
    }
  }

  void SoftwareRendererImp::rasterize_image(float x0, float y0,
    float x1, float y1,
    Texture& tex) {
    // Task 4: 
    // Implement image rasterization (you may want to call fill_sample here)

  }

  // resolve samples to render target
  void SoftwareRendererImp::resolve(void) {

    // Task 2: 
    // Implement supersampling
    // You may also need to modify other functions marked with "Task 2".
    return;

  }


} // namespace CS248
