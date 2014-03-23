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
#include "Matrix4x4.h"
using FasTC::Matrix4x4;
using FasTC::Vec4f;

#include "Pixel.h"
using FasTC::Pixel;
using FasTC::YCoCgPixel;

#include "Image.h"
#include "ImageFile.h"
#include "StopWatch.h"
#include "CompressedImage.h"

#include "BPTCCompressor.h"

#include "SLIC.h"
#include "Partition.h"
#include "VPTree.h"

static Matrix4x4<float> ComputeCovarianceMatrix(
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
  
  Matrix4x4<float> covMatrix;

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
      if(fabs(fabs(v * dir) - v.Length()) > 1e-7) {
        collinear = false;
        break;
      }
    }

    if(collinear) {
      axis = dir;
      return 0;
    }
  }

  return ComputeCovarianceMatrix(pts).PowerMethod(axis);
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

static YCoCgPixel Vec4fToPixel(const Vec4f &v) {
  YCoCgPixel p;
  p.Y() = CastChannel(Clamp(0.0f, 255.0f, v[0]));
  p.Co() = CastChannel(Clamp(0.0f, 255.0f, v[1]));
  p.Cg() = CastChannel(Clamp(0.0f, 255.0f, v[2]));
  p.A() = CastChannel(Clamp(0.0f, 255.0f, v[3]));
  return p;
}

static void PrintPixel(const char *label, const Pixel &p) {
  fprintf(stderr, "%s: <%d, %d, %d, %d>\n", label, p.R(), p.G(), p.B(), p.A());
}

class Region {
  YCoCgPixel m_Endpoints[2];
  std::vector<Pixel>::const_iterator m_PixelItr;
  std::vector<Pixel> m_Pixels;
  std::vector<uint8> m_Interp;
  std::vector<uint8> m_LumaInterp;

  double CompressVecs(std::vector<Vec4f> pts) {

//    for(const auto &p : pts) {
//      fprintf(stderr, "%s: <%.2f, %.2f, %.2f, %.2f>\n", "Pt", p.X(), p.Y(), p.Z(), p.W());
//    }

    std::vector<float> firstChannel;
    firstChannel.reserve(pts.size());
    for(auto &p : pts) {
      firstChannel.push_back(p[0]);
      p[0] = 255.0f;
    }

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

    // Luma...
    float fmin = *min_element(std::begin(firstChannel), std::end(firstChannel));
    float fmax = *max_element(std::begin(firstChannel), std::end(firstChannel));
    m_Endpoints[0].Y() = static_cast<int16>(fmin);
    m_Endpoints[1].Y() = static_cast<int16>(fmax);

    m_LumaInterp.reserve(m_Pixels.size());
    m_Interp.reserve(m_Pixels.size());
    for(uint32 i = 0; i < m_Pixels.size(); i++) {
      
      auto p = pts[i];
      float f = firstChannel[i];

      if(b == a) {
        m_Interp.push_back(0);
      } else {
        float d = (p - centroid).Dot(axis);
        float nd = (d - a) / (b - a);
        assert(0.0f <= nd && nd <= 1.0f);
        m_Interp.push_back(FloatToChannel(nd));
      }

      if(fmax == fmin) {
        m_LumaInterp.push_back(0);
      } else {
        float d = (f - fmin) / (fmax - fmin);
        assert(0.0f <= d && d <= 1.0f);
        m_LumaInterp.push_back(FloatToChannel(d));
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

    assert(opaque);
//    if(opaque) {

      m_Endpoints[0][3] = m_Endpoints[1][3] = 255;

      for(const auto &p : m_Pixels) {
        YCoCgPixel trans = YCoCgPixel(p);
        Vec4f v;
        v[0] = static_cast<float>(static_cast<int16>(trans.Y()));
        v[1] = static_cast<float>(static_cast<int16>(trans.Co()));
        v[2] = static_cast<float>(static_cast<int16>(trans.Cg()));
        v[3] = 255.0f;
        rgbaVecs.push_back(v);
      }
//    } else {
//
//      for(const auto &p : m_Pixels) {
//        Vec4f v;
//        v[0] = static_cast<float>(p.R());
//        v[1] = static_cast<float>(p.G());
//        v[2] = static_cast<float>(p.B());
//        v[3] = static_cast<float>(p.A());
//        rgbaVecs.push_back(v);
//      }
//    }

    return CompressVecs(rgbaVecs);
  }

  void Reconstruct() {
    m_Pixels.clear();
    m_Pixels.reserve(m_Interp.size());

    for(uint32 i = 0; i < m_Interp.size(); i++) {
      const auto &v = m_Interp[i];
      const auto &y = m_LumaInterp[i];

      float yd = static_cast<float>(y) / 255.0f;
      float vd = static_cast<float>(v) / 255.0f;

      YCoCgPixel p = m_Endpoints[0] * (1 - vd) + m_Endpoints[1] * vd;
      p.Y() = m_Endpoints[0].Y() * (1 - yd) + m_Endpoints[1].Y() * yd;
//      PrintPixel("Reconstructed", p);
      m_Pixels.push_back(p.ToRGBA());
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

typedef VpTree<Partition<4, 4>, Partition<4, 4>::Distance> VpTree4x4;
struct SelectionInfo {
  VpTree4x4 &tree;
  const int *labels;
  const uint32 width;
  const uint32 height;

  SelectionInfo(VpTree4x4 &_t, const int *l, uint32 w, uint32 h)
    : tree(_t)
    , labels(l)
    , width(w)
    , height(h)
  { }
};

static uint32 kTwoPartitionModes = 
  static_cast<uint32>(BPTCC::eBlockMode_One) |
  static_cast<uint32>(BPTCC::eBlockMode_Three) |
  static_cast<uint32>(BPTCC::eBlockMode_Seven);
static uint32 kThreePartitionModes =
  static_cast<uint32>(BPTCC::eBlockMode_Zero) |
  static_cast<uint32>(BPTCC::eBlockMode_Two);
#if 0
static uint32 kAlphaModes =
  static_cast<uint32>(BPTCC::eBlockMode_Four) |
  static_cast<uint32>(BPTCC::eBlockMode_Five) |
  static_cast<uint32>(BPTCC::eBlockMode_Six)  |
  static_cast<uint32>(BPTCC::eBlockMode_Seven);
#endif

template<const unsigned N, const unsigned M>
BPTCC::ShapeSelection ChosePresegmentedShape(
  uint32 x, uint32 y, const uint32 pixels[16], const void *userData
) {
  const SelectionInfo &info = *(reinterpret_cast<const SelectionInfo *>(userData));

  // Construct a partition...
  Partition<N, M> part;

  static const uint32 kMaxLabelsPerShape = 6;
  int32 map[kMaxLabelsPerShape];
  memset(map, 0xFF, sizeof(map));

  bool opaque = true;

  uint32 idx = 0;
  for(uint32 i = x; i < x+N; i++)
  for(uint32 j = y; j < y+M; j++) {
    int label = info.labels[j*info.width + i];

    // Has this label been seen already?
    uint32 l = 0;
    for(; l < kMaxLabelsPerShape; l++) {
      if(map[l] == label) {
        break;
      } else if(map[l] < 0) {
        map[l] = label;
        break;
      }
    }

    assert(l < kMaxLabelsPerShape);
    part[idx] = l;

    if(((pixels[idx] >> 24) & 0xFF) < 250) {
      opaque = false;
    }
  }

  // Is this a single partition?
  BPTCC::ShapeSelection result;
  if(map[1] < 0) {
    // Turn off two and three shape modes.
    result.m_SelectedModes &=
      ~(kThreePartitionModes | kTwoPartitionModes);
  } else {

    std::vector<Partition<N, M> > closestParts;
    info.tree.search(part, 1, &closestParts, NULL);

    const Partition<N, M> &closest = closestParts[0];
    uint8 maxPart = 0;
    for(uint32 i = 0; i < N*M; i++) {
      maxPart = std::max(maxPart, closest[i]);
    }

    if(maxPart < 2) {
      result.m_TwoShapeIndex = closest.GetIndex();
      // Turn off three shape modes
      result.m_SelectedModes &= ~(kThreePartitionModes);
    } else {
      result.m_ThreeShapeIndex = closest.GetIndex();
      // Turn off two shape modes
      result.m_SelectedModes &= ~(kTwoPartitionModes);      
    }
  }

  // If opaque, turn off modes 4 and five...
  if(opaque) {
    result.m_SelectedModes &=
      ~(static_cast<uint32>(BPTCC::eBlockMode_Four) |
        static_cast<uint32>(BPTCC::eBlockMode_Five));
  }

  return result;
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

  std::vector<Partition<4, 4> > partitions;
  EnumerateBPTC(partitions);
  std::cout << partitions.size() << " 4x4 BPTC partitions" << std::endl;

  VpTree<Partition<4, 4>, Partition<4, 4>::Distance> vptree;
  vptree.create(partitions);

  // Just to test, find the partition close to half 0 half 1..
  Partition<4, 4> test;
  for(uint32 i = 0; i < 16; i++) {
    if(i < 8) {
      test[i] = 0;
    } else {
      test[i] = 1;
    }
  }

  vector<Partition<4, 4> > closest;
  vptree.search(test, 1, &closest, NULL);
  std::cout << closest[0].GetIndex() << std::endl;

  BPTCC::CompressionSettings settings;
  settings.m_NumSimulatedAnnealingSteps = 0;
  settings.m_ShapeSelectionFn = ChosePresegmentedShape<4, 4>;

  SelectionInfo info(vptree, labels, kWidth, kHeight);
  settings.m_ShapeSelectionUserData = &info;

  uint8 *outBuf = new uint8[kWidth * kHeight];
  FasTC::CompressionJob cj(
     FasTC::eCompressionFormat_BPTC,
     reinterpret_cast<const uint8 *>(pixels),
     outBuf,
     static_cast<uint32>(kWidth),
     static_cast<uint32>(kHeight));

  StopWatch sw;
  sw.Start();
  BPTCC::Compress(cj, settings);
  sw.Stop();
  std::cout << "Compression time: " << sw.TimeInMilliseconds() << "ms" << std::endl;

  CompressedImage ci(kWidth, kHeight, FasTC::eCompressionFormat_BPTC, outBuf);
  FasTC::Image<> outImg(kWidth, kHeight, pixels);

  std::cout << "PSNR: " << outImg.ComputePSNR(&ci) << "db" << std::endl;

  ImageFile outImgFile("out.png", eFileFormat_PNG, outImg);
  outImgFile.Write();

  delete [] labels;
  delete [] rawPixels;
  delete [] pixels;
  return 0;
}
