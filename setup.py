from distutils.core import setup

setup(name="camb_utils",
      version="0.0.1dev",
      description='Utilities releated to reading camb files',
      long_description=''' ''',
      packages=['camb_utils'],
      package_dir={'camb_utils': 'camb_utils'},
      include_package_data=True,
      package_data={'camb_utils': ['example_data/*']}
      

      )
