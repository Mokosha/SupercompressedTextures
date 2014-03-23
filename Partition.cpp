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

#include "Partition.h"

#include <cassert>
#include <cstdlib>
#include <cstring>

// Partition selection functions as specified in
// C.2.21
static inline uint32 hash52(uint32 p) {
  p ^= p >> 15;  p -= p << 17;  p += p << 7; p += p << 4;
  p ^= p >>  5;  p += p << 16;  p ^= p >> 7; p ^= p >> 3;
  p ^= p <<  6;  p ^= p >> 17;
  return p;
}

static int32 SelectPartition(int32 seed, int32 x, int32 y, int32 z,
                             int32 partitionCount, int32 smallBlock) {
  if(smallBlock) {
    x <<= 1;
    y <<= 1;
    z <<= 1;
  }

  seed += (partitionCount-1) * 1024;

  uint32 rnum = hash52(seed);
  uint8 seed1  =  rnum        & 0xF;
  uint8 seed2  = (rnum >>  4) & 0xF;
  uint8 seed3  = (rnum >>  8) & 0xF;
  uint8 seed4  = (rnum >> 12) & 0xF;
  uint8 seed5  = (rnum >> 16) & 0xF;
  uint8 seed6  = (rnum >> 20) & 0xF;
  uint8 seed7  = (rnum >> 24) & 0xF;
  uint8 seed8  = (rnum >> 28) & 0xF;
  uint8 seed9  = (rnum >> 18) & 0xF;
  uint8 seed10 = (rnum >> 22) & 0xF;
  uint8 seed11 = (rnum >> 26) & 0xF;
  uint8 seed12 = ((rnum >> 30) | (rnum << 2)) & 0xF;

  seed1 *= seed1;     seed2 *= seed2;
  seed3 *= seed3;     seed4 *= seed4;
  seed5 *= seed5;     seed6 *= seed6;
  seed7 *= seed7;     seed8 *= seed8;
  seed9 *= seed9;     seed10 *= seed10;
  seed11 *= seed11;   seed12 *= seed12;

  int32 sh1, sh2, sh3;
  if(seed & 1) {
    sh1 = (seed & 2)? 4 : 5;
    sh2 = (partitionCount == 3)? 6 : 5;
  } else {
    sh1 = (partitionCount == 3)? 6 : 5;
    sh2 = (seed & 2)? 4 : 5;
  }
  sh3 = (seed & 0x10) ? sh1 : sh2;

  seed1 >>= sh1; seed2  >>= sh2; seed3  >>= sh1; seed4  >>= sh2;
  seed5 >>= sh1; seed6  >>= sh2; seed7  >>= sh1; seed8  >>= sh2;
  seed9 >>= sh3; seed10 >>= sh3; seed11 >>= sh3; seed12 >>= sh3;

  int32 a = seed1*x + seed2*y + seed11*z + (rnum >> 14);
  int32 b = seed3*x + seed4*y + seed12*z + (rnum >> 10);
  int32 c = seed5*x + seed6*y + seed9 *z + (rnum >>  6);
  int32 d = seed7*x + seed8*y + seed10*z + (rnum >>  2);

  a &= 0x3F; b &= 0x3F; c &= 0x3F; d &= 0x3F;

  if( partitionCount < 4 ) d = 0;
  if( partitionCount < 3 ) c = 0;

  if( a >= b && a >= c && a >= d ) return 0;
  else if( b >= c && b >= d ) return 1;
  else if( c >= d ) return 2;
  return 3;
}

template<uint32 N, uint32 M>
void EnumerateASTC(std::vector<Partition<N, M> > &results) {
  uint8 parts[144]; // Largest possible block size
  const uint32 smallBlock = (N*M) < 31;

  results.clear();
  
  uint32 kMaxPartitionIndex = (1 << 12) - 1;
  for(uint32 partIdx = 0; partIdx <= kMaxPartitionIndex; partIdx++) {
    memset(parts, 0, sizeof(parts));
    for(uint32 j = 0; j < M; j++)
    for(uint32 i = 0; i < N; i++) {
      uint8 part = SelectPartition(partIdx, i, j, 0, partIdx & 0x3, smallBlock);
      assert(part < 4);
      parts[j*N + i] = part;
    }

    // Remap so that the parts are in increasing order
    int32 map[4];
    memset(map, 0xFF, sizeof(map));

    int32 lastPart = 0;
    for(uint32 j = 0; j < M; j++)
    for(uint32 i = 0; i < N; i++) {
      uint8 &part = parts[j*N + i];
      if(map[part] < 0) {
        map[part] = lastPart++;
      }
      part = map[part];
    }
    assert(lastPart <= 4);

    Partition<N, M> partition(partIdx);
    for(uint32 i = 0; i < N*M; i++) {
      partition[i] = parts[i];
    }

    // Make sure that the partition doesn't already exist
    bool exists = false;
    for(uint32 i = 0; i < results.size(); i++) {
      if(results[i] == partition) {
        exists = true;
        break;
      }
    }

    if(!exists) {
      results.push_back(partition);
    }
  }
}

// Anything else will cause a link error...
template void EnumerateASTC<4, 4>(std::vector<Partition<4, 4> > &);
template void EnumerateASTC<5, 4>(std::vector<Partition<5, 4> > &);
template void EnumerateASTC<5, 5>(std::vector<Partition<5, 5> > &);
template void EnumerateASTC<6, 5>(std::vector<Partition<6, 5> > &);
template void EnumerateASTC<6, 6>(std::vector<Partition<6, 6> > &);
template void EnumerateASTC<8, 5>(std::vector<Partition<8, 5> > &);
template void EnumerateASTC<8, 6>(std::vector<Partition<8, 6> > &);
template void EnumerateASTC<8, 8>(std::vector<Partition<8, 8> > &);
template void EnumerateASTC<10, 5>(std::vector<Partition<10, 5> > &);
template void EnumerateASTC<10, 6>(std::vector<Partition<10, 6> > &);
template void EnumerateASTC<10, 8>(std::vector<Partition<10, 8> > &);
template void EnumerateASTC<10, 10>(std::vector<Partition<10, 10> > &);
template void EnumerateASTC<12, 10>(std::vector<Partition<12, 10> > &);
template void EnumerateASTC<12, 12>(std::vector<Partition<12, 12> > &);

#include "Shapes.h"

void EnumerateBPTC(std::vector<Partition<4, 4> > &results) {
  uint8 parts[16]; // Largest possible block size
  results.clear();
  
  uint32 kMaxShapeIndex = 63;
  for(uint32 nSubsets = 2; nSubsets <= 3; nSubsets++)
  for(uint32 shapeIdx = 0; shapeIdx <= kMaxShapeIndex; shapeIdx++) {
    memset(parts, 0, sizeof(parts));
    for(uint32 i = 0; i < 16; i++) {
      uint8 part = BPTCC::GetSubsetForIndex(i, shapeIdx, nSubsets);
      assert(part < 3);
      parts[i] = part;
    }

    // Remap so that the parts are in increasing order
    int32 map[4];
    memset(map, 0xFF, sizeof(map));

    int32 lastPart = 0;
    for(uint32 i = 0; i < 16; i++) {
      uint8 &part = parts[i];
      if(map[part] < 0) {
        map[part] = lastPart++;
      }
      part = map[part];
    }
    assert(lastPart <= 4);

    Partition<4, 4> partition(shapeIdx);
    for(uint32 i = 0; i < 16; i++) {
      partition[i] = parts[i];
    }

    // Make sure that the partition doesn't already exist
    bool exists = false;
    for(uint32 i = 0; i < results.size(); i++) {
      if(results[i] == partition) {
        exists = true;
        break;
      }
    }

    if(!exists) {
      results.push_back(partition);
    }
  }  
}

