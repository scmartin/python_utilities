import setuptools

with open("README.md", "r") as fh:
	    long_description = fh.read()

setuptools.setup(
    name="autocorr",
    version="0.0.1",
    author="Seth Martin",
    description="Contains code for calculating autocorrelation \
                 of time series data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

