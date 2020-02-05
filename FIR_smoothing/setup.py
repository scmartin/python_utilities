import setuptools

with open("README.md", "r") as fh:
	    long_description = fh.read()

setuptools.setup(
    name="FIR_smoothing",
    version="0.0.1",
    author="Seth Martin",
    description="Contains module with functions for FIR_smoothing fine scale \
    profiles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

