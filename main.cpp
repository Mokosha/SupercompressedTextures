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

#include "Image.h"
#include "ImageFile.h"

#ifdef _MSC_VER
int _tmain(int argc, _TCHAR* argv[]) {
#else
int main(int argc, char **argv) {
#endif

  if(argc != 3) {
    fprintf(stderr, "Usage: compare <img1> <img2>\n");
    return 1;
  }

  ImageFile img1f (argv[1]);
  if(!img1f.Load()) {
    fprintf(stderr, "Error loading file: %s\n", argv[1]);
    return 1;
  }

  ImageFile img2f (argv[2]);
  if(!img2f.Load()) {
    fprintf(stderr, "Error loading file: %s\n", argv[2]);
    return 1;
  }

  FasTC::Image<> img1(*img1f.GetImage());
  FasTC::Image<> img2(*img2f.GetImage());

  double PSNR = img1.ComputePSNR(&img2);
  if(PSNR > 0.0) {
    fprintf(stdout, "PSNR: %.3f\n", PSNR);
  }
  else {
    fprintf(stderr, "Error computing PSNR\n");
  }

  double SSIM = img1.ComputeSSIM(&img2);
  if(SSIM > 0.0) {
    fprintf(stdout, "SSIM: %.9f\n", SSIM);
  } else {
    fprintf(stderr, "Error computing MSSIM\n");
  }

  return 0;
}
