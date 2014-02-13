#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include <unordered_map>
#ifdef _MSC_VER
#  include <SDKDDKVer.h>
#  include <Windows.h>
#endif

#include "Vector4.h"
#include "MatrixSquare.h"
using FasTC::MatrixSquare;
typedef FasTC::Vector4<float> Vec4f;

#include "Pixel.h"
using FasTC::Pixel;

#include "Image.h"
#include "ImageFile.h"

#include "SLIC.h"

static MatrixSquare<float, 4> ComputeCovarianceMatrix(
  const std::vector<Vec4f> &points
) {
  
  Vec4f avg;
  for(int i = 0; i < 4; i++)
    avg[i] = float(0.0f);

  for(const auto &p : points) {
    avg += p;
  }
  avg /= float(points.size());

  std::vector<Vec4f > toPts;
  toPts.reserve(points.size());

  for(const auto &p : points) {
    toPts.push_back(p - avg);
  }
  
  MatrixSquare<float, 4> covMatrix;

  // Compute covariance.
  for(uint32 i = 0; i < 4; i++) {
    for(uint32 j = 0; j <= i; j++) {

      float sum(0.0);
      for(uint32 k = 0; k < points.size(); k++) {
        sum += toPts[k][i] * toPts[k][j];
      }

      covMatrix(i, j) = sum / float(4 - 1);
      covMatrix(j, i) = covMatrix(i, j);
    }
  }

  return covMatrix;
}

struct CompareVecs {
  bool operator()(const Vec4f &v1, const Vec4f &v2) {
    for(uint32 i = 0; i < 4; i++) {
      if(v1[i] != v2[i]) {
        return v1[i] < v2[i];
      }
    }
    return false;
  }
};

static uint32 GetPrincipalAxis(const std::vector<Vec4f> &pts,
                               Vec4f &axis) {

  // Generate a list of unique points...
  CompareVecs cv;
  std::set<Vec4f, CompareVecs> upts(cv);
  for(const auto &p : pts) {
    upts.insert(p);
  }

  assert(upts.size() > 0);

  if(upts.size() == 1) {
    for(int i = 0; i < 4; i++)
      axis[i] = float(0.0);
    return 0;

  // Collinear?
  } else {
    typedef typename std::set<Vec4f >::const_iterator Itr;
    Itr itr = upts.begin();
    assert(itr != upts.end());

    Vec4f upt0 = *itr;
    itr++;
    assert(itr != upts.end());

    Vec4f dir (*itr - upt0);
    dir.Normalize();
    itr++;

    bool collinear = true;
    for(; itr != upts.end(); itr++) {
      Vec4f v = (*itr - upt0);
      if(fabs(fabs(v.Dot(dir)) - v.Length()) > 1e-7) {
        collinear = false;
        break;
      }
    }

    if(collinear) {
      axis = dir;
      return 0;
    }
  }

  MatrixSquare<float, 4> covMatrix = ComputeCovarianceMatrix(pts); 
  uint32 iters = covMatrix.PowerMethod(axis);

#if 0
  if(NULL != eigTwo) {
    if(eigOne != 0.0) {
      RGBAMatrix reduced = covMatrix - eigOne * RGBAMatrix(
        axis.c[0] * axis.c[0], axis.c[0] * axis.c[1], axis.c[0] * axis.c[2], axis.c[0] * axis.c[3], 
        axis.c[1] * axis.c[0], axis.c[1] * axis.c[1], axis.c[1] * axis.c[2], axis.c[1] * axis.c[3],
        axis.c[2] * axis.c[0], axis.c[2] * axis.c[1], axis.c[2] * axis.c[2], axis.c[2] * axis.c[3],
        axis.c[3] * axis.c[0], axis.c[3] * axis.c[1], axis.c[3] * axis.c[2], axis.c[3] * axis.c[3]
      );
      
      bool allZero = true;
      for(uint32 i = 0; i < 16; i++) {
        if(fabs(reduced[i]) > 0.0005) {
          allZero = false;
        }
      }

      if(allZero) {
        *eigTwo = 0.0;
      }
      else {
        RGBADir dummyDir;
        iters += PowerIteration(reduced, dummyDir, *eigTwo);
      }
    }
    else {
      *eigTwo = 0.0;
    }
  }
#endif

  return iters;
}

template<typename T>
static inline T Clamp(const T &a, const T &b, const T &v) {
  return std::max(a, std::min(v, b));
}

static uint8 CastChannel(const float &f) {
  return static_cast<uint8>(f + 0.5f);
}

static uint8 FloatToChannel(const float &f) {
  return CastChannel(255.0f * f);
}

static Pixel Vec4fToPixel(const Vec4f &v) {
  Pixel p;
  p.R() = CastChannel(Clamp(0.0f, 255.0f, v[0]));
  p.G() = CastChannel(Clamp(0.0f, 255.0f, v[1]));
  p.B() = CastChannel(Clamp(0.0f, 255.0f, v[2]));
  p.A() = CastChannel(Clamp(0.0f, 255.0f, v[3]));
  return p;
}

static void PrintPixel(const char *label, const Pixel &p) {
  fprintf(stderr, "%s: <%d, %d, %d, %d>\n", label, p.R(), p.G(), p.B(), p.A());
}

class Region {
  Pixel m_Endpoints[2];
  std::vector<Pixel>::const_iterator m_PixelItr;
  std::vector<Pixel> m_Pixels;
  std::vector<uint8> m_InterpolationValues;

  double CompressVecs(const std::vector<Vec4f> &pts) {

//    for(const auto &p : pts) {
//      fprintf(stderr, "%s: <%.2f, %.2f, %.2f, %.2f>\n", "Pt", p.X(), p.Y(), p.Z(), p.W());
//    }

    Vec4f centroid = Vec4f(0, 0, 0, 0);
    for(const auto &p : pts) {
      centroid += p;
    }
    centroid /= pts.size();
//    fprintf(stderr, "Centroid: <%.2f, %.2f, %.2f, %.2f>\n", centroid.X(), centroid.Y(), centroid.Z(), centroid.W());

    float a = FLT_MAX, b = -FLT_MAX;
    Vec4f axis;
    if(GetPrincipalAxis(pts, axis) <= 0) {
      m_Endpoints[0] = m_Endpoints[1] = Vec4fToPixel(centroid);
      a = b = 0.0f;
    } else {
      axis.Normalize();

      for(const auto &pt : pts) {
        float d = (pt - centroid).Dot(axis);
        a = std::min(d, a);
        b = std::max(d, b);
      }

      m_Endpoints[0] = Vec4fToPixel(centroid + (axis * a));
      m_Endpoints[1] = Vec4fToPixel(centroid + (axis * b));
    }

    m_InterpolationValues.reserve(m_Pixels.size());
    for(const auto &pt : pts) {
      if(b == a) {
        m_InterpolationValues.push_back(0);
      } else {
        float d = (pt - centroid).Dot(axis);
        float nd = (d - a) / (b - a);
        assert(0.0f <= nd && nd <= 1.0f);
      
        m_InterpolationValues.push_back(FloatToChannel(nd));
      }
    }

    return b - a;
  }

public:
  Region() { }
  explicit Region(const std::vector<Pixel> &pixels)
    : m_Pixels(pixels) { }

  uint32 NumPixels() const { return m_Pixels.size(); }
  void AddPixel(const Pixel &p) {
    m_Pixels.push_back(p);
  }

  void ResetPixelItr() { m_PixelItr = m_Pixels.begin(); }
  Pixel GetNextPixel() {
    assert(m_PixelItr != m_Pixels.end());
    const Pixel &ret = *m_PixelItr;
    m_PixelItr++;
    return ret;
  }

  // Populate m_Endpoints and m_InterpolationValues
  double Compress() {
    bool opaque = true;
    for(const auto &p : m_Pixels) {
      opaque = opaque && (p.A() >= 255);
    }

    std::vector<Vec4f> rgbaVecs;
    rgbaVecs.reserve(m_Pixels.size());

    if(opaque) {

      m_Endpoints[0].A() = m_Endpoints[1].A() = 255;

      for(const auto &p : m_Pixels) {
        Vec4f v;
        v[0] = static_cast<float>(p.R());
        v[1] = static_cast<float>(p.G());
        v[2] = static_cast<float>(p.B());
        v[3] = 255.0f;
        rgbaVecs.push_back(v);
      }
    } else {

      for(const auto &p : m_Pixels) {
        Vec4f v;
        v[0] = static_cast<float>(p.R());
        v[1] = static_cast<float>(p.G());
        v[2] = static_cast<float>(p.B());
        v[3] = static_cast<float>(p.A());
        rgbaVecs.push_back(v);
      }
    }

    return CompressVecs(rgbaVecs);
  }

  void Reconstruct() {
    m_Pixels.clear();
    m_Pixels.reserve(m_InterpolationValues.size());
    for(const auto &v : m_InterpolationValues) {
      float nd = static_cast<float>(v) / 255.0f;
      Pixel p = m_Endpoints[0] * (1 - nd) + m_Endpoints[1] * nd;
//      PrintPixel("Reconstructed", p);
      m_Pixels.push_back(p);
    }

    ResetPixelItr();
  }
};

void CollectPixels(const uint32 kWidth, const uint32 kHeight, 
                   const Pixel *pixels, const int *labels,
                   std::unordered_map<uint32, Region> &result) {
  result.clear();
  for(uint32 j = 0; j < kHeight; j++) {
    for(uint32 i = 0; i < kWidth; i++) {
      uint32 idx = j*kWidth+i;
      uint32 label = static_cast<uint32>(labels[idx]);
      Pixel p = pixels[idx];

      if(result.count(label) == 0) {
        Region r;
        r.AddPixel(p);
        std::pair<uint32, Region> newRegion (label, r);
        result.insert(newRegion);
      } else {
        result[label].AddPixel(p);
      }
    }
  }
}

#ifdef _MSC_VER
int _tmain(int argc, _TCHAR* argv[]) {
#else
int main(int argc, char **argv) {
#endif

  if(argc != 2 && argc != 3) {
    fprintf(stderr, "Usage: sc <img1>\n");
    return 1;
  }

  int spSize = 5;
  if(argc == 3) {
    sscanf(argv[2], "%d", &spSize);
  }

  ImageFile imgFile (argv[1]);
  if(!imgFile.Load()) {
    fprintf(stderr, "Error loading file: %s\n", argv[1]);
    return 1;
  }

  FasTC::Image<> *img = imgFile.GetImage();

  const int kWidth = img->GetWidth();
  const int kHeight = img->GetHeight();
  const int nPixels = kWidth * kHeight;
  const uint32 pixelBufSz = nPixels * sizeof(FasTC::Pixel);

  FasTC::Pixel *pixels = new FasTC::Pixel[pixelBufSz];
  memcpy(pixels, img->GetPixels(), pixelBufSz);

  uint32 *rawPixels = new uint32[kWidth * kHeight];

  for(int i = 0; i < nPixels; i++) {
    // Pixels are stored as little endian ARGB, so we want ABGR
    pixels[i].Shuffle(0x6C); // 01 10 11 00
    rawPixels[i] = pixels[i].Pack();
  }

  int *labels = new int[nPixels];
  int numLabels;

  SLIC slic;
  slic.PerformSLICO_ForGivenStepSize(
    rawPixels,
	kWidth,
    kHeight,
    labels,
    numLabels,
	spSize, 1.0);

  std::unordered_map<uint32, Region> regions;
  CollectPixels(kWidth, kHeight, pixels, labels, regions);
  std::cout << "Num regions: " << regions.size() << std::endl;

  for(auto &r : regions) {
    r.second.Compress();
    r.second.Reconstruct();
  }

  for(int i = 0; i < nPixels; i++) {
    pixels[i] = regions[labels[i]].GetNextPixel();
    pixels[i].Shuffle(0x6C);
  }

  FasTC::Image<> outImg(kWidth, kHeight, pixels);
  ImageFile outImgFile("out.png", eFileFormat_PNG, outImg);
  outImgFile.Write();

  delete [] labels;
  delete [] rawPixels;
  delete [] pixels;
  return 0;
}
