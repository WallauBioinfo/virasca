import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="virseqimprover",
    version="1.2",
    author="Haoqiu Song",
    author_email="haoqiu@vt.edu",
    description="An integrated pipeline for error-correction, extension, and annotation of viral scaffolds.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/haoqiusong/VirSeqImprover",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
    ),
)
