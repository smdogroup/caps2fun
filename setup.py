"""
Setup file for caps2fun package
run "pip install -e ." in this dir to install package
"""

from setuptools import setup

setup(
    name="caps_to_funtofem",
    version="0.1",
    author="Sean Engelstad",
    author_email="sengelstad312@gatech.edu",
    description="interface between ESP/CAPS, Funtofem, TACS, and FUN3D for Aerothermoelastic analysis of Aerospace Structures",
    packages=["pycaps2fun"],
    python_requires=">=3.6",
)