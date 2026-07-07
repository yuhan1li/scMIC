import setuptools


def get_readme():
    with open("README.md", "rt", encoding="utf-8") as handle:
        return handle.read()


setuptools.setup(
    name="motmic",
    version="0.1.0",
    author="yuhan1li",
    description="multi-organ optimal transport for metastasis-initiating cell discovery",
    long_description=get_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yuhan1li/scMIC",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy>=1.21",
        "scipy>=1.7",
        "pandas>=1.3",
        "scikit-learn>=1.0",
        "matplotlib>=3.5",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 3 - Alpha",
    ],
    python_requires=">=3.8",
)

