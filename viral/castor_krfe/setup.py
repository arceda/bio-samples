from setuptools import setup, find_packages
from castor_krfe import __version__ as _version

INSTALL_REQUIRES = []

with open("requirements.txt", "r") as fh:
    for line in fh:
        INSTALL_REQUIRES.append(line.rstrip())

setup(
    name='castor_krfe',
    version=_version,
    description='',
    author='Dylan',
    author_email='lebatteux.dylan@courrier.uqam.ca',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'castor_krfe = Main:main'
            ]
        },
    install_requires=INSTALL_REQUIRES
)
