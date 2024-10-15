"""Setup and dependencies"""

from setuptools import setup, find_packages

setup(
    name="datafunks",
    version="0.1.0",
    author="MF",
    description="These modules are collections of scripts I've used over the years to perform basic analysis tasks in the database.",
    packages=find_packages(),
    install_requires=[
        # Project dependencies
        "numpy",
        "pandas",
        "matplotlib",
        "shapely",
        "geopandas",
        "holoviews",
        "datashader",
        "colorcet",
        "scipy",
        "lifelines",
        "statsmodels"
    ],
    python_requires="<3.11",  # Should be 3.10 to work with astropathdb
    entry_points={
        "console_scripts": [
            # Add command line scripts here if needed
        ],
    },
)

# For ipynbs- instead of restaring kernel each time
# %load_ext autoreload
# %autoreload 2
