#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#ifdef _MSC_VER
#  include <SDKDDKVer.h>
#  include <Windows.h>
#endif

#include "Pixel.h"
using FasTC::Pixel;

#include "Image.h"
#include "ImageFile.h"

#include "SLIC.h"

class Region {
  Pixel m_Endpoints[2];
  std::vector<Pixel> m_Pixels;
  std::vector<uint8> m_InterpolationValues;

public:
  Region() { }
  explicit Region(const std::vector<Pixel> &pixels)
    : m_Pixels(pixels) { }

  uint32 NumPixels() const { return m_Pixels.size(); }
  void AddPixel(const Pixel &p) {
    m_Pixels.push_back(p);
  }

  // Populate m_Endpoints and m_InterpolationValues
  void Compress() {
  }
};

void CollectPixels(const uint32 kWidth, const uint32 kHeight, 
                   const Pixel *pixels, const int *labels,
                   std::unordered_map<uint32, Region> &result) {
  result.clear();
  for(uint32 j = 0; j < kWidth; j++) {
    for(uint32 i = 0; i < kHeight; i++) {
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

  if(argc != 2) {
    fprintf(stderr, "Usage: sc <img1>\n");
    return 1;
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
	6, 1.0);

  std::unordered_map<uint32, Region> regions;
  CollectPixels(kWidth, kHeight, pixels, labels, regions);
  std::cout << "Num regions: " << regions.size() << std::endl;

  uint32 np = 0;
  for(const auto &rp : regions) {
    uint32 num = rp.second.NumPixels();
    np += num;
  }
  std::cout << "Total pixels in regions: " << np << std::endl;

  for(int i = 0; i < nPixels; i++) {
    pixels[i].Unpack(rawPixels[i]);
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
