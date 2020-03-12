#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Extra credit: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  MipLevel& mip = tex.mipmap[0];
  int width = mip.width;
  int height = mip.height;
  int x = u * width;
  if (x == width) {
    x = width - 1;
  }
  int y = v * height;
  if (y == height) {
    y = height - 1;
  }

  int start = 4 * (y * mip.width + x);
  Color color;
  uint8_to_float(&color.r, &(mip.texels[start]));
  return color;

}

static void getIndex(float v, int length, int& floor, int& ceil)
{
  if (v < 0.5) {
    floor = 0;
    ceil = 0;
  }
  else if (v >(length - 0.5)) {
    floor = length - 1;
    ceil = length - 1;
  }
  else {
    floor = static_cast<int>(v - 0.5);
    ceil = floor + 1;
  }
}
static Color lerp(float s, Color v1, Color v2) {
  return v1 + s * (v2 - v1);
}

Color Sampler2DImp::getColor(Texture& tex, int x, int y) {
  MipLevel& mip = tex.mipmap[0];
  int width = mip.width;
  int height = mip.height;
  int index = 4 * (y * mip.width + x);
  Color color;
  uint8_to_float(&color.r, &(mip.texels[index]));
  return color;
}
Color Sampler2DImp::sample_bilinear(Texture& tex,
                                    float u, float v, 
                                    int level) {
  MipLevel& mip = tex.mipmap[0];
  int width = mip.width;
  int height = mip.height;
  float x = u * width;
  float y = v * height;
  int floorX;
  int ceilX;
  getIndex(x, width, floorX, ceilX);
  int floorY;
  int ceilY;
  getIndex(y, height, floorY, ceilY);
  
  Color topLeft= getColor(tex, floorX, floorY);
  Color topRight = getColor(tex, ceilX, floorY);
  Color bottomLeft = getColor(tex, floorX, ceilY);
  Color bottomRight = getColor(tex, ceilX, ceilY);

  float sx = x - (floorX + 0.5);
  Color colorUp = lerp(sx, topLeft, topRight);
  Color colorDown = lerp(sx, bottomLeft, bottomRight);
  float sy = y - (floorY + 0.5);
  return lerp(sy, colorUp, colorDown);
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Extra credit: Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248
