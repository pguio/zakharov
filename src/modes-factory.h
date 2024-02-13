/**************************************************************************
 *
 * $Id: modes-factory.h,v 1.6 2011/03/26 07:47:28 patrick Exp $
 *
 * Copyright (c) 2004-2011 Patrick Guio <patrick.guio@gmail.com>
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


#ifndef SPECIES_FACTORY_H
#define SPECIES_FACTORY_H

#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <classexception.h>


typedef std::string defaultIDKeyType;

namespace factory {
  void dummySolitons();
  void dummyEfluct();
  void dummyEgauss();
  void dummyThermalFluct();
  void dummyDensityFluct();
  void dummyDensityParab();
}

template <class ManufacturedType, typename ClassIDKey=defaultIDKeyType>
class ModesFactory {
public:
  typedef ManufacturedType* BasePtr;
private:

  typedef BasePtr (*BaseCreateFn)(int , char **);
  typedef std::map<ClassIDKey, BaseCreateFn> FnRegistry;
  FnRegistry registry;
  ModesFactory()
  {}
  ModesFactory(const ModesFactory&); // Not implemented
  ModesFactory &operator=(const ModesFactory&); // Not implemented
public:
  typedef ClassIDKey IDKeyType;
  typedef typename FnRegistry::const_iterator const_iterator;

  static ModesFactory &instance() {
    static ModesFactory<ManufacturedType, ClassIDKey> bf;
    return bf;
  }

  void RegCreateFn(const ClassIDKey & id, BaseCreateFn fn) {
    registry[id] = fn;
  }

  BasePtr create(const ClassIDKey &className,
                 int nargs, char *args[]) const {
    BasePtr theObject(0);
    const_iterator regEntry = registry.find(className);
    if (regEntry != registry.end()) {
      theObject = regEntry->second(nargs,args);
    } else {
      throw ClassException("ModesFactory", "Error: unknown ClassIDKey: " + className);
    }
    return theObject;
  }
  const_iterator begin() const {
    return registry.begin();
  }
  const_iterator end  () const {
    return registry.end();
  }
  ClassIDKey getKey(const_iterator &i) const {
    return i->first;
  }

  void content() const {
    using std::cout;
    using std::endl;

    using std::left;
    using std::hex;
    using std::right;
    using std::setw;

    using std::string;

    const_iterator i = registry.begin();
    cout << "Factory content: " << endl;
    cout << string(30, '=') << endl;
    cout << setw(20) << left << "Name"
         << setw(10) << right << "Address" << endl;
    cout << string(30, '=') << endl;
    for (const_iterator i = registry.begin(); i != registry.end(); i++) {
      cout << setw(20) << left << i->first
           << setw(10) << right << hex << int(i->second)  << endl;
    }
  }
  void init() const {
    factory::dummySolitons();
    factory::dummyEfluct();
    factory::dummyEgauss();
    factory::dummyThermalFluct();
    factory::dummyDensityFluct();
    factory::dummyDensityParab();
  }
};



template <class AncestorType,
         class ManufacturedType,
         typename ClassIDKey=defaultIDKeyType>
class RegisterInFactory {
private:
  typedef ModesFactory<AncestorType, ClassIDKey> Factory;
  typedef typename Factory::BasePtr BasePtr;
public:
  static BasePtr CreateInstance(int nargs, char *args[]) {
    return BasePtr(new ManufacturedType(nargs, args));
  }

  RegisterInFactory(const ClassIDKey &id) {
    Factory::instance().RegCreateFn(id, CreateInstance);
  }
};


#endif
