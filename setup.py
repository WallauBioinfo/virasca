from setuptools import setup, find_packages

setup(
    name="virasca",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "Click",
        "snakemake==7.32.4",
        "pulp<2.8",
        "pandas",
        "pyyaml"
    ],
    entry_points={
        "console_scripts": [
            "virasca = virasca.cli:cli",
        ],
    },
)
