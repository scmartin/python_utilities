import setuptools

with open("README.md", "r") as fh:
	    long_description = fh.read()

setuptools.setup(
    name="lammps",
    version="0.0.1",
    author="Seth Martin",
    description="Contains modules for loading lammps logs and trajectories",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

