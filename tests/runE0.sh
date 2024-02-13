#!/bin/sh

#
# $Id: runE0.sh,v 1.2 2011/03/26 08:30:21 patrick Exp $
#
# Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
# All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2.  of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

dir=../hdf/beamE0
filename=beam

case $HOSTTYPE in
	i386 | *linux) 
    timelogger="/usr/bin/time -v -o run.log --append"
		;;
  solaris)
    timelogger="/usr/bin/time -p"
		;;
	alpha)
    timelogger="/usr/bin/time"
		;;
	*)
    echo 'A timelogger has to be defined'
		exit 0
esac

execdir=`pwd`
execname=$execdir/Zakharov
conffile=$execdir/../share/beam.conf

if [ ! -d $dir ] ; then
	install -d $dir
fi

cd $dir

if [ -f run.log ]; then
	rm run.log
fi

E0=(1e-4 7.5e-5 5e-5 2.5e-5 1e-5 7.5e-6 5e-6 2.5e-6 1e-6 0)
ext=(1e-4 7.5e-5 5e-5 2.5e-5 1e-5 7.5e-6 5e-6 2.5e-6 1e-6 0)
iter="0 1 2 3 4 5 6 7 8 9"

for i in ${iter} ; do
	command="$execname -o $filename.hdf -i $conffile E0=${E0[i]} avSkw2nbk=2"
#	echo $timelogger $command 
	echo $timelogger $command > run.cmd
 	($timelogger $command > run.log) > run.cmd 2>&1
 	mv $filename.hdf $filename${ext[i]}.hdf
# 	echo mv $filename.hdf $filename${ext[i]}.hdf
 	mv run.cmd $filename${ext[i]}.cmd
 	mv run.log $filename${ext[i]}.log
	gzip -f9 $filename${ext[i]}.log	
done
cp $conffile .
