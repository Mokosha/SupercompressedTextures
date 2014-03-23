/* FasTC
 * Copyright (c) 2014 University of North Carolina at Chapel Hill.
 * All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and non-profit purposes, without
 * fee, and without a written agreement is hereby granted, provided that the
 * above copyright notice, this paragraph, and the following four paragraphs
 * appear in all copies.
 *
 * Permission to incorporate this software into commercial products may be
 * obtained by contacting the authors or the Office of Technology Development
 * at the University of North Carolina at Chapel Hill <otd@unc.edu>.
 *
 * This software program and documentation are copyrighted by the University of
 * North Carolina at Chapel Hill. The software program and documentation are
 * supplied "as is," without any accompanying services from the University of
 * North Carolina at Chapel Hill or the authors. The University of North
 * Carolina at Chapel Hill and the authors do not warrant that the operation of
 * the program will be uninterrupted or error-free. The end-user understands
 * that the program was developed for research purposes and is advised not to
 * rely exclusively on the program for any reason.
 *
 * IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL OR THE
 * AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
 * THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF NORTH CAROLINA
 * AT CHAPEL HILL OR THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 * THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL AND THE AUTHORS SPECIFICALLY
 * DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY 
 * STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON
 * AN "AS IS" BASIS, AND THE UNIVERSITY  OF NORTH CAROLINA AT CHAPEL HILL AND
 * THE AUTHORS HAVE NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
 * ENHANCEMENTS, OR MODIFICATIONS.
 *
 * Please send all BUG REPORTS to <pavel@cs.unc.edu>.
 *
 * The authors may be contacted via:
 *
 * Pavel Krajcevski
 * Dept of Computer Science
 * 201 S Columbia St
 * Frederick P. Brooks, Jr. Computer Science Bldg
 * Chapel Hill, NC 27599-3175
 * USA
 * 
 * <http://gamma.cs.unc.edu/FasTC/>
 */

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "VPTree.h"
#include "TexCompTypes.h"
#include "StopWatch.h"

static double Hamming(const uint32 &a, const uint32 &b) {
  uint32 n = a ^ b;
  uint32 c;
  for(c = 0; n; c++) {
    n &= n-1;
  }
  return static_cast<double>(c);
}

static const uint32 kNumVals = 4096;

int main() {
  srand(time(NULL));

  std::vector<uint32> vals;
  vals.reserve(kNumVals);
  for(uint32 i = 0; i < kNumVals; i++) {
    vals.push_back(rand());
  }

  StopWatch stopwatch;

  const uint32 target = 0x0000beef;
  uint32 bestFive[5] = {~target, ~target, ~target, ~target, ~target};
  std::cout << "Linear: " << std::endl;
  stopwatch.Start();
  for(unsigned int i = 0; i < vals.size(); i++) {
    for(uint32 j = 0; j < 5; j++) {
      if(Hamming(target, vals[i]) < Hamming(target, bestFive[j])) {
        for(uint32 k = 4; k > j; k--) {
          bestFive[k] = bestFive[k - 1];
        }
        bestFive[j] = vals[i];
        break;
      }
    }
  }
  stopwatch.Stop();
  std::cout << "Time: (" << stopwatch.TimeInMilliseconds() << " ms)" << std::endl;
  for(unsigned int i = 0; i < 5; i++) {
    std::cout << std::dec << i << " (" << Hamming(bestFive[i], target) << "): 0x" << std::hex << bestFive[i] << std::endl;
  }
  std::cout << std::endl;

  VpTree<uint32, Hamming> vptree;
  vptree.create(vals);

  std::cout << "VPTree: " << std::endl;
  stopwatch.Reset();
  stopwatch.Start();
  std::vector<uint32> results;
  for(uint32 i = 0; i < 100; i++) {
    vptree.search(target, 1, &results, NULL);
  }
  stopwatch.Stop();
  std::cout << "Time: (" << (stopwatch.TimeInMilliseconds() / 100.0) << " ms)" << std::endl;
  for(unsigned int i = 0; i < results.size(); i++) {
    std::cout << std::dec << i << " (" << Hamming(results[i], target) << "): 0x" << std::hex << results[i] << std::endl;
  }
  std::cout << std::endl;
}
