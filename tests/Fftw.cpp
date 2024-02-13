/**************************************************************************
 *
 * $Id: Fftw.cpp,v 1.13 2011/03/26 08:30:21 patrick Exp $
 *
 * Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
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

#include <cstdlib>
#include <iostream>
using namespace std;

#include <fourier.h>
using namespace blitz;

int main(int nargs, char *args[])
{
  try {
    int n=64;

    if (nargs == 2) n = atoi(args[1]);

    fourier::ODFT1D<fourier::complex,fourier::complex> ft1d;
    Array<complex<double>, 1> A(n), B(n), C(n);

    ft1d.resize(n,-1);
    cout << "ft1d = " << ft1d << endl;
    A=1.0;
    cout << "A=" << A << endl;
    ft1d.direct(A, B);
    cout << "B=" << B << endl;
    ft1d.inverse(B, C);
    cout << "C=" << C << endl;

    cout << "B(0)=" << B(0) << " sum(A)=" << sum(A) << endl;

    ft1d.resize(n,1);
    cout << "ft1d = " << ft1d << endl;
    A=1.0;
    cout << "A=" << A << endl;
    ft1d.direct(A, B);
    cout << "B=" << B << endl;
    ft1d.inverse(B, C);
    cout << "C=" << C << endl;

    cout << "B(0)=" << B(0) << " sum(A)=" << sum(A) << endl;

    return EXIT_SUCCESS;

  } catch (int status) {
    return status;
  } catch (ClassException& c) {
    cerr << c.what() << endl;
    return EXIT_FAILURE;
  }

}
