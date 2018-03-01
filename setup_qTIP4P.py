# Use this script to compile the Cython module by running
#
#         python setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import os, sys

try:
  os.remove("MMTK_qTIP4P.so")
except:
  print 'dont need to delete'
  
compile_args = []
include_dirs = ['/home/kpbishop/old_Dev_GPU/Downloads/khinsen-mmtk-35865b3a97f5/Include']

from Scientific import N
try:
  num_package = N.package
except AttributeError:
  num_package = "Numeric"
if num_package == "NumPy":
  compile_args.append("-DNUMPY=1")
  import numpy.distutils.misc_util
  include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())
  
setup (name = "qTIP4PforceField",
       version = "1.0",
       description = "q-TIP4P forcefield term for MMTK",
       
       py_modules = ['qTIP4pFF'],
       ext_modules = [Extension('MMTK_qTIP4p',
                                ['MMTK_qTIP4p.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs)],
       cmdclass = {'build_ext': build_ext}
     )
