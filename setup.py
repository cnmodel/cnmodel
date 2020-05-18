from setuptools import setup, find_packages
import setuptools
import os, shutil


import os

path = os.path.join(os.path.dirname(__file__), 'cnmodel')
version = None
for line in open(os.path.join(path, '__init__.py'), 'r').readlines():
    if line.startswith('__version__'):
        version = line.partition('=')[2].strip('"\' \n')
        break
if version is None:
    raise Exception("Could not read __version__ from cnmodel/__init__.py")

setup(name='cnmodel',
      version=version,
      description='A biophysically realistic model of cochlear nucleus neurons and other neurons',
      url='http://github.com/pbmanis/cnmodel',
      author='Paul B. Manis and Luke Campagnola',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['cnmodel*']),
      install_requires=['matplotlib>=3.0', 'numpy>=1.1',
          ],
      entry_points={
          'console_scripts': [
               'CNtoy_model=examples.toy_model:main',
               'CNtest_cells=examples.test_cells:main',
               'CNunittests=test:main',
               ],

      },
      zip_safe=False,
      data_files=[('mechs', ['x86_64/*'])],  # includes the current compiled mechanisms
#      cmdclass={'makeneuron': 'Build_Nmodl'},
      classifiers = [
             "Programming Language :: Python :: 3.7+",
             "Development Status ::  Beta",
             "Environment :: Console",
             "Intended Audience :: Neuroscientists, computational",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Software Development :: Tools :: Python Modules",
             "Topic :: Computational Modeling :: Neuroscience",
             ],
      )
      