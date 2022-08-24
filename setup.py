"""
Setup file for caps2fun package
run "pip install -e ." in this dir to install package
"""

from setuptools import setup

setup(
    name="caps2fun",
    version="1.0",
    author="Sean Engelstad",
    author_email="sengelstad312@gatech.edu",
    description="interface between ESP/CAPS, Funtofem, TACS, and FUN3D for Aerothermoelastic analysis of Aerospace Structures",
    packages=[
        "capsManager", "tacsManager", "fun3dManager", "funtofemManager",
        "caps2tacs", "caps2fun", "caps2funtofem", "utils"
        ],
    package_dir = {'': 'src'},
    python_requires=">=3.6",
)