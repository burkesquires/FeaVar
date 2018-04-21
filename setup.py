from setuptools import setup, find_packages

setup(name='ombre',

      version='0.1',

      url='https://github.com/burkesquires/ombre',

      license='MIT',

      author='R. Burke Squires',

      author_email='burkesquires@gmail.com',

      description='A python package to compute clusters of sequence feature variant types (SFVTs) \
      based upon user-selected subsequence',

      packages=find_packages(exclude=['tests']),

      long_description=open('README.md').read(),

      zip_safe=False,

      setup_requires=['nose>=1.0'],

      test_suite='nose.collector',

      install_requires=['biopython',
                        'pandas']
      )