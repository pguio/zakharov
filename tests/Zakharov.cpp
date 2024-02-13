/**************************************************************************
 *
 * $Id: Zakharov.cpp,v 1.27 2019/05/10 16:55:11 patrick Exp $
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

#include <new>
using namespace std;

#include <zakharov.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

// Default format parameters
const int defaultPrecision = 7;
const std::ios::fmtflags defaultFloatfield = std::ios::scientific;


void my_new_handler()
{
  cerr << "Out of memory" << endl;
  abort();
}

#define PARSE(FUN)             \
if (parser.FUN()) {            \
	solver.FUN();                \
	return 0;                    \
}                              \


int main(int nargs, char *args[])
{
  try {

    std::set_new_handler(&my_new_handler);
    cout.precision(defaultPrecision);
    cout.setf(defaultFloatfield, std::ios::floatfield);
#if 0
    printOsStatus(cout);
#endif

#if defined(HAVE_MPI)
    MPI_Init(&nargs, &args);
#endif

    parser::Parser parser(nargs, args);
    parser.registerProgram(args[0]);
    parser.registerPackage(PACKAGE, VERSION, ZAK_COPYRIGHT);

    ZakharovSolver solver(nargs, args);

    PARSE(parseHelp)
    PARSE(parseVersion)
    PARSE(parseTemplate)

    solver.initialise();

    cout << solver << endl << endl;

    solver.solve();

#if defined(HAVE_MPI)
    MPI_Finalize();
#endif

    return 0;

  } catch (const ClassException& c) {
    cerr << c.what() << endl;
    return !0;
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
    return !0;

  } catch (const int status) {
    cerr << "status=" << status << endl;
    return status;
  }
}


