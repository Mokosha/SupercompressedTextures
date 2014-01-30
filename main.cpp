#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#ifdef _MSC_VER
#  include <SDKDDKVer.h>
#  include <Windows.h>
#endif

#include "Pixel.h"
using FasTC::Pixel;

#include "Image.h"
#include "ImageFile.h"

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

  slic.DrawContoursAroundSegments(
	rawPixels,
	labels,
	kWidth,
	kHeight,
	0xFF000000);

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
