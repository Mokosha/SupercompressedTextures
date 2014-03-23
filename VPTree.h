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

#ifndef _VPTREE_H__
#define _VPTREE_H__ 

// A VP-Tree implementation, by Steve Hanov. (steve.hanov@gmail.com)
// Released to the Public Domain
// Based on "Data Structures and Algorithms for Nearest Neighbor Search" by Peter N. Yianilos
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>

template<typename T, double (*distance)( const T&, const T& )>
class VpTree
{
 private:
  template<class _T>
  class ReservablePQueue : public std::priority_queue<_T> {
   public:
    typedef typename std::priority_queue<T>::size_type size_type;
    ReservablePQueue(size_type capacity = 0) { reserve(capacity); };
    void reserve(size_type capacity) { this->c.reserve(capacity); } 
    size_type capacity() const { return this->c.capacity(); } 
    void clear() { this->c.clear(); }
  };
  
 public:
 VpTree() : _root(0), _heap(20) {}

  ~VpTree() {
    delete _root;
  }

  void create( const std::vector<T> &items ) {
    delete _root;
    _items = items;
    _root = buildFromPoints(0, items.size());
  }

  void search( const T& target, int k, std::vector<T> *results, 
               std::vector<double> *distances) 
  {
    _heap.clear();
    _tau = std::numeric_limits<double>::max();
    search( _root, target, k );

    results->clear();
    results->reserve(k);
    if(distances) {
      distances->clear();
      distances->reserve(k);
    }

    while( !_heap.empty() ) {
      results->push_back( _items[_heap.top().index] );
      if(distances)
        distances->push_back( _heap.top().dist );
      _heap.pop();
    }

    std::reverse( results->begin(), results->end() );
    if(distances)
      std::reverse( distances->begin(), distances->end() );
  }

 private:
  std::vector<T> _items;
  double _tau;

    struct Node 
    {
      int index;
      double threshold;
      Node* left;
      Node* right;

    Node() :
      index(0), threshold(0.), left(0), right(0) {}

      ~Node() {
        delete left;
        delete right;
      }
    }* _root;

    struct HeapItem {
    HeapItem( int index, double dist) :
      index(index), dist(dist) {}
      int index;
      double dist;
      bool operator<( const HeapItem& o ) const {
        return dist < o.dist;   
      }
    };
    ReservablePQueue<HeapItem> _heap;

    struct DistanceComparator
    {
      const T& item;
    DistanceComparator( const T& item ) : item(item) {}
      bool operator()(const T& a, const T& b) {
        return distance( item, a ) < distance( item, b );
      }
    };

    Node* buildFromPoints( int lower, int upper )
    {
      if ( upper == lower ) {
        return NULL;
      }

      Node* node = new Node();
      node->index = lower;

      if ( upper - lower > 1 ) {

        // choose an arbitrary point and move it to the start
        int i = (int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;
        std::swap( _items[lower], _items[i] );

        int median = ( upper + lower ) / 2;

        // partitian around the median distance
        std::nth_element( 
                         _items.begin() + lower + 1, 
                         _items.begin() + median,
                         _items.begin() + upper,
                         DistanceComparator( _items[lower] ));

        // what was the median?
        node->threshold = distance( _items[lower], _items[median] );

        node->index = lower;
        node->left = buildFromPoints( lower + 1, median );
        node->right = buildFromPoints( median, upper );
      }

      return node;
    }

    void search( Node* node, const T& target, int k )
    {
      if ( node == NULL ) return;

      double dist = distance( _items[node->index], target );
      //printf("dist=%g tau=%gn", dist, _tau );

      if ( dist < _tau ) {
        if ( _heap.size() == k ) _heap.pop();
        _heap.push( HeapItem(node->index, dist) );
        if ( _heap.size() == k ) _tau = _heap.top().dist;
      }

      if ( node->left == NULL && node->right == NULL ) {
        return;
      }

      if ( dist < node->threshold ) {
        if ( dist - _tau <= node->threshold ) {
          search( node->left, target, k );
        }

        if ( dist + _tau >= node->threshold ) {
          search( node->right, target, k );
        }

      } else {
        if ( dist + _tau >= node->threshold ) {
          search( node->right, target, k );
        }

        if ( dist - _tau <= node->threshold ) {
          search( node->left, target, k );
        }
      }
    }
};

#endif // _VPTREE_H__
