#!/bin/sh

case `hostname` in

phalanx*)
	./configure --with-plplot=/usr/include/plplot,/usr/lib64 --with-fftw3=/usr/include,/usr/lib64 --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-blitz=$HOME/src/blitz/,$HOME/src/blitz/lib/.libs --disable-cxx-flags-preset CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-std=c++0x -xSSE4.1 -O2 -restrict -vec-report1 -no-prec-div -no-ansi-alias" AR=xiar
	;;

alun*|dare*) 
#	./configure --with-fftw3=/usr/include,/usr/lib64 --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-blitz=$HOME --enable-cxx-flags-preset --enable-optimize CXX="g++"
	./configure --with-fftw3 --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-blitz=$HOME --disable-cxx-flags-preset --enable-optimize CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-xSSSE3 -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias" AR=xiar
	;;
			
kryten*|hawke*) 
#	./configure --with-fftw=$HOME --with-hdf4=$HOME --with-blitz=$HOME --enable-cxx-flags-preset --enable-optimize CXX="g++"
	./configure --with-fftw=$HOME --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-blitz=$HOME --enable-cxx-flags-preset --enable-optimize CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-xSSSE3 -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias" AR=xiar
	;;
			
eukalyptus*)
#./configure --with-fftw3=/net/pl/fftw --with-hdf4=/net/pl/hdf --with-blitz=/net/pl/blitz --disable-cxx-flags-preset CXX=g++ CXXFLAGS="-g -O2" 
./configure --with-fftw3=/net/pl/fftw --with-hdf4=/net/pl/hdf --with-blitz=/net/pl/blitz --enable-optimize CXX=g++  
;;

fimm*)
./configure  --with-blitz=$HOME --with-fftw3=/local/fftw --with-hdf4=/local/HDF  --enable-optimize CXX=g++
#./configure  --with-blitz=$HOME --with-fftw3=/local/fftw --with-hdf4=$HOME/src/hdf --enable-optimize CXX=g++
;;

login*)
        if test "$LOGNAME" = "ucappgu"; then
        ./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --disable-cxx-flags-preset --enable-optimize CXX="icpc" CXXFLAGS="-O3 -ip -static -no-prec-div -xP -axP"
        elif test "$LOGNAME" = "patrickg"; then
#       ./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --enable-cxx-flags-preset --enable-optimize CXX="g++" CXXFLAGS="-O2"
#       ./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-optimize CXX="pathCC" CXXFLAGS="-Ofast"
        ./configure --with-fftw3=/site/fftw3.intel --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-optimize --enable-threadsafe-blitz  CXX="icc -Kc++" CXXFLAGS="-xW -axW -O3 -ip -static -no-prec-div -ansi-alias"
#       ./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-optimize CXX="pgCC" CXXFLAGS="-fast"
  fi
        ;;

jazz5*|jazz6*)
./configure CONFIG_SITE=config/sites/f32_x86_64_g++_optim_jazz5
;;

jazz2*|jazz3*)
./configure CONFIG_SITE=config/sites/f31_x86_64_g++_optim_jazz3
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local/blitz1 --with-plplot=/usr/include/plplot,/usr/lib --disable-cxx-flags-preset CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-std=c++0x -xSSE4.2 -O3 -restrict -vec-report=2 -no-prec-div -no-ansi-alias" AR=xiar
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local/blitz2 --with-plplot=/usr/include/plplot,/usr/lib --disable-cxx-flags-preset CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-std=c++0x -xSSE4.2 -O3 -restrict -vec-report=2 -no-prec-div -no-ansi-alias" AR=xiar
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --with-plplot=/usr/include/plplot,/usr/lib --disable-cxx-flags-preset CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-std=c++0x -xSSE4.2 -O3 -ipo -restrict -vec-report=2 -no-prec-div -no-ansi-alias" AR=xiar
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --with-plplot=/usr/include/plplot,/usr/lib --disable-cxx-flags-preset CXX=icpc --enable-threadsafe-blitz CXXFLAGS="-xSSSE3 -O3 -ipo -restrict -vec-report1 -no-prec-div -no-ansi-alias" AR=xiar
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --with-plplot=/usr/include/plplot,/usr/lib --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-std=c++0x -g -fp-stack-check -DBZ_DEBUG"
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-plplot=/usr/include/plplot,/usr/lib --with-blitz --enable-optimize CXX=g++
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz --disable-cxx-flags-preset CXX=g++ CXXFLAGS="-Wall -pedantic -O3 -funroll-loops -fstrict-aliasing -fomit-frame-pointer -ffast-math -DBZ_DEBUG"
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz --disable-cxx-flags-preset CXX=g++ CXXFLAGS="-O2"
#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz --disable-cxx-flags-preset CXX=g++ CXXFLAGS=""
;;

jazz*|fyspc-rp42*|theta*)
#./configure --with-blitz=/usr/local --with-fftw=/usr/local/fftw --with-hdf4=/usr/local/hdf CXX=icc
./configure --with-fftw3 --with-hdf4=/usr/include/hdf,/usr/lib/hdf --enable-gnuplot --disable-cxx-flags-preset CXX=g++ CXXFLAGS="-O2"
#./configure --with-fftw3 --with-blitz=/usr/local --with-hdf4=/usr/include/hdf,/usr/lib/hdf --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-xN -O3 -ip -no-prec-div -ansi-alias"
#./configure gcc --with-blitz=/usr/local --with-fftw=/usr/local/fftw --with-hdf4=/usr/local/hdf CXX=g++
#./configure --with-blitz=/usr/local --with-fftw=/usr/local/fftw --with-hdf4=/usr/local/hdf --enable-debug CXX=g++
;;

phact*|hyades*|ayil*|barakish*|menkar*|markab*|papsukal*|thuban*)
./configure --with-blitz=$HOME --with-fftw3=$HOME --with-hdf4=$HOME/src/hdf --enable-dependency-tracking  --disable-cxx-flags-preset --enable-optimize CC=cc CXX=cxx CXXFLAGS="-std ansi -D__USE_STD_IOSTREAM -DBZ_ENABLE_XOPEN_SOURCE -D_OSF_SOURCE -ieee -model ansi -accept restrict_keyword -nousing_std -no_implicit_include" CXX_OPTIMIZE_FLAGS="-fast -inline speed -nocleanup"
#./configure CXX=cxx --with-blitz=$HOME --enable-dxml --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --enable-debug
#./configure CXX=cxx --with-blitz=$HOME --enable-dxml --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --enable-optimize
;;

tre*)
	./configure CXX=xlC --enable-64bit --with-blitz=$HOME --with-fftw=/usr/local/fftw-2.1.5-64 --with-hdf4=/usr/local/HDF4-64 --enable-optimize
	;;

gridur*|embla*|balder*)
	./configure CC=cc --with-cxx=SGI64 --with-blitz=$HOME --with-fftw=/usr/local/lib --with-hdf4=/usr/local/hdf --enable-optimize
	;;
	
magnum*|pico*)
	./configure --with-cxx=aCC --with-blitz=$HOME --with-fftw=/site/fftw --with-hdf4=/site/hdf --enable-dependency-tracking --enable-optimize
	;;

*titan*)
#	./configure --with-fftw=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --enable-optimize CXX=g++
#	./configure --with-fftw=$HOME/include,$HOME/lib/intel --with-blitz=$HOME/include,$HOME/lib/intel --with-hdf4=$HOME/include,$HOME/lib/intel --enable-optimize CXX=icpc 
#	./configure --with-fftw=$HOME/include,$HOME/lib/intel --with-blitz=$HOME/include,$HOME/lib/intel --with-hdf4=$HOME/include,$HOME/lib/intel --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-xP -O3 -ip -no-prec-div -ansi-alias"
#	./configure --with-fftw=$HOME/include,$HOME/lib/portland --with-blitz=$HOME/include,$HOME/lib/portland --with-hdf4=$HOME/include,$HOME/lib/portland --enable-optimize CXX=pgCC
	./configure --with-fftw=$HOME/include,$HOME/lib/portland --with-blitz=$HOME/include,$HOME/lib/portland --with-hdf4=$HOME/include,$HOME/lib/portland --disable-cxx-flags-preset CXX=pgCC CXXFLAGS="-fastsse"
#	./configure --with-fftw=$HOME/include,$HOME/lib/portland --with-blitz=$HOME/include,$HOME/lib/portland --with-hdf4=$HOME/include,$HOME/lib/portland --disable-cxx-flags-preset CXX=pgCC CXXFLAGS=" -O4 -Mnoframe -Mnodepchk -Minline"
	;;

snowstorm*)
	./configure --with-cxx=icc --with-blitz=$HOME --with-fftw=$HOME --with-hdf=$HOME --enable-dependency-tracking --enable-optimize
	;;

nana*)
	#./configure --with-cxx=aCC --with-blitz=$HOME --with-fftw=/usr/local/numerics/fftw --with-hdf4=/usr/local/hdf/hdf --enable-dependency-tracking --enable-optimize
	./configure --with-cxx=aCC --with-blitz=$HOME --with-mlib=/opt/mlib --with-hdf4=/usr/local/hdf/hdf --enable-dependency-tracking --enable-optimize
	#./configure --with-cxx=aCC --with-blitz=$HOME --with-fftw=/usr/local/numerics/fftw --with-hdf4=/usr/local/hdf/hdf --enable-dependency-tracking --enable-optimize
	;;


*) echo No default configuration for machine `hostname`

esac
