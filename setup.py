from setuptools import setup 
  
# reading long description from file 
with open('DESCRIPTION.txt') as file: 
    long_description = file.read() 
# specify requirements of your package here 
REQUIREMENTS = ['requests'] 
  
# some more details 
CLASSIFIERS = [ 
    'Development Status :: 4 - Beta', 
    'Intended Audience :: Developers', 
    'Topic :: Scientific/Engineering :: Chemistry', 
    'Programming Language :: Python', 
    'Programming Language :: Python :: 2', 
    'Programming Language :: Python :: 2.6', 
    'Programming Language :: Python :: 2.7', 
    'Programming Language :: Python :: 3', 
    'Programming Language :: Python :: 3.3', 
    'Programming Language :: Python :: 3.4', 
    'Programming Language :: Python :: 3.5', 
    ] 
  
# calling the setup function  
setup(name='cg2at', 
      version='v0.2', 
      description='Fragment based conversion of the CG representation into atomistic.', 
      long_description=long_description, 
      url='https://github.com/owenvickery/cg2at', 
      author='Owen Vickery', 
      author_email='owen.vickery@warwick.ac.uk', 
      license='GPLv3', 
      classifiers=CLASSIFIERS, 
      install_requires=REQUIREMENTS, 
      keywords='cg2at'
      ) 

