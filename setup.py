from setuptools import setup, find_packages

setup(
    name="kinfraglib",
    version="2.0.0",
    description="Kinase-focused fragment library",
    url="https://github.com/volkamerlab/kinfraglib",
    author="Volkamer Lab",
    author_email="volkamer@cs.uni-saarland.de",
    license="MIT",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
)
