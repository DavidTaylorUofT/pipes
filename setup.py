'''setup for cython module defined in allthethings.pyx 
Compile with command
	CXX=g++ python setup.py build_ext -i

Use in python with command:
	from allthethings import [stuff] (so far [stuff] can be: * or PyNetwork or PyBC_opt_dh or PyMystery_BC)

Runs on my macbook air, my desktop Mac, and on a linux VM. I'm doing something kind of sketch with the lapack wrapper to get around segfaults that happen if I call the libcla.a library like in the pure C++ code in Build/. by kind of sketch I mean...I copied lapack.c and added it as a dependency, and wrapped the header.
BASICALLY ALL OF THIS IS BLACK MAGIC... 
...so if it doesn't compile for you, I'm terribly sorry for the hassle and I will try to help...but no guarantees :/
!!!!!PS if it's yelling about w convert 64 to 32 not recognized by gcc4.9, try running these two commands in the shell before compiling
(or just type "source cython_flags.sh" in your terminal):

    $ export ARCHFLAGS=""           
    $ export CFLAGS="-arch i386 -arch x86_64" 
'''
'''This file is part of Pipes.

    Pipes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pipes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pipes.  If not, see <http://www.gnu.org/licenses/>.
'''


import sys
import os

##!!!!! below: you should specify path to anaconda distribution and gcc compiler.
#I uncomment either the line below (1.1), (1.2), or (1.3)

# (1.1) on my desktop:
#sys.path.append('/Users/lieba/anaconda/lib/python2.7/site-packages')
# (1.2) on my macbook air
#sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages')
#sys.path.append('/usr/local/Cellar/gcc49/4.9.2_1/lib/gcc/4.9/gcc/x86_64-apple-darwin12.6.0/4.9.2/include-fixed')
#(1.3)if on linux VM (replace [USERNAME] with your username. Verify this location with $ which python2.7)
sys.path.append('/home/xin/anaconda2/lib/python2.7/site-packages')


##!!!! below: you should specify compiler including path.
# I uncomment either the lines below (2.1),(2.2), or (2.3)
#(2.1) on my desktop
#os.environ["CC"] = "gcc-4.9"
#os.environ["CXX"] = "g++-4.9"
#(2.2) on my macbook air
#os.environ["CC"] = "/usr/local/bin/gcc-4.9"
#os.environ["CXX"] = "/usr/local/bin/gcc-4.9"
#(2.3) on a linux VM
os.environ["CC"] = "/usr/bin/gcc-4.9"
os.environ["CXX"] = "/usr/bin/g++-4.9"

from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy
print 'start'
setup(ext_modules=cythonize(Extension(
    "allthethings",                                   # the extesion name
    sources=["allthethings.pyx", "channel.cpp", "setupandrun.cpp", "file_output.cc", "network.cpp", "levmar.cpp", "mp_mat.cpp",
             "str_double.cpp", "mp_mat_double.cpp", "libcla.c","setup.py"],  # the Cython source and additional C++ source files
    # sources=["allthethings.pyx", "setupandrun.cpp", "file_output.cc","network.cpp", "levmar.cpp","mp_mat.cpp","str_double.cpp", "mp_mat_double.cpp"], # the Cython source and additional C++ source files
    # libraries=["lapack","cblas", "qd", "fftw3","m"],      #libraries to link
    # against (I'm not sure if all are needed...but some are...)
    # libraries to link against (I'm not sure if all are needed...but some
    # are...)
    libraries=["lapack", "cblas", "fftw3", "m"],
    language="c++",                         # generate and compile C++ code

    # on orinoco
    # extra_link_args=['-DUSEOMP'],
    extra_link_args=['-fopenmp'],
    # extra_compile_args=['-DUSEOMP'],
    extra_compile_args=['-fopenmp'],
    # so it can find, e.g. numpy/arrayobject.h
    # on orinoco
    include_dirs=[numpy.get_include(),"/home/xin", "/usr/include"]   #so it can find, e.g. numpy/arrayobject.h
    # on macbook Air
    #include_dirs=[numpy.get_include(), "/Users/anna", "/usr/local/include"]
)))
print 'end'
