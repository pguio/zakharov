#!/bin/sh

#
# $Id: egauss.sh,v 1.2 2011/03/26 08:30:21 patrick Exp $
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

prefix=/mn/sothi/u1/patrickg/astro-cluster/research/codes/zakharov

#L=2400
#L=1200
L=600

filename=egauss
dir=hdf/${filename}${L}
execname=tests/Zakharov
conffile=share/${filename}${L}.conf

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


if [ ! -d $dir ] ; then
	install -d $dir
fi

if [ -f run.log ]; then
	rm run.log
fi

for Erms in 1e-2 1e-3 1e-4 1e-5 1e-6 0e-6 ; do
	for nrms in 1e-4 1e-5 1e-6 1e-7 1e-8 0e-8 ; do
		if [ ! -f $dir/$filename$Erms$nrms.hdf ] ; then
			command="$execname Erms=$Erms nrms=$nrms -o $filename.hdf -i $conffile"
			echo $timelogger $command > run.cmd
 			($timelogger $command > run.log) > run.cmd 2>&1
 			mv $filename.hdf $dir/$filename$Erms$nrms.hdf
 			mv run.cmd $dir/$filename$Erms$nrms.cmd
 			mv run.log $dir/$filename$Erms$nrms.log
			gzip -f9 $dir/$filename$Erms$nrms.log	
		fi
	done
done
cp $conffile $dir
