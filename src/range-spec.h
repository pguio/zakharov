/**************************************************************************
 *
 * $Id: range-spec.h,v 1.3 2011/11/08 11:47:48 patrick Exp $
 *
 * Copyright (c) 2010-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef RANGE_SPEC_H
#define RANGE_SPEC_H

#if !defined(BLITZ2)
#include <blitz/tinyvec-et.h>
#else
#include <blitz/tinyvec2.h>
#endif
#include <blitz/range.h>
#include <iostream>


template <class P_numtype>
class RangeSpec {
public:

  typedef blitz::TinyVector<P_numtype,3> Triplet;

  friend std::ostream& operator<<(std::ostream &os, const RangeSpec & r) {
    if (r.size() == 1)
      return os << r.start();
    else
      return os << r.start() << ":" << r.stride() << ":" << r.end();
  }

  RangeSpec(P_numtype val=0) : data(val)
  {}
  RangeSpec(Triplet val) : data(val)
  {}
  RangeSpec(P_numtype start, P_numtype stride, P_numtype end) :
    data(start, stride, end)
  {}

  P_numtype start()  const {
    return data(0);
  }
  P_numtype stride() const {
    return data(1);
  }
  P_numtype end()    const {
    return data(2);
  }

  void start(P_numtype _x) {
    data(0) = _x;
  }
  void stride(P_numtype _x) {
    data(1) = _x;
  }
  void end(P_numtype _x) {
    data(2) = _x;
  }

  P_numtype pos(P_numtype val)  const {
    return ( val - this->start() ) / this->stride();
  }

  int size() const {
    return static_cast<int>((this->end()-this->start())/this->stride() + 1);
  }

  bool isWithin(P_numtype val) const {
    return val >= this->start() && val <= this->end();
  }

  bool operator!=(P_numtype val) const {
    return blitz::all(data != val);
  }
  bool operator!=(const RangeSpec<P_numtype> &val) const {
    return blitz::all(data != val);
  }
  bool operator==(P_numtype val) const {
    return blitz::all(data == val);
  }
  bool operator==(const RangeSpec<P_numtype> &val) const {
    return blitz::all(data == val);
  }

  Triplet data;

};

inline
blitz::Range toRange(const RangeSpec<int> &spec)
{
  return blitz::Range(spec.start(), spec.end(), spec.stride());
}


#endif
