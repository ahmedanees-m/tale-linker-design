from setuptools import setup, find_packages

setup(
    name="tale-linker-design",
    version="0.1.0",
    description="Geometric framework for TALE-fusion catalytic domain linker design",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Anees Ahmed",
    author_email="ahmedaneesm@gmail.com",
    url="https://github.com/ahmedanees-m/tale-linker-design",
    license="MIT",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.24",
        "scipy>=1.10",
        "pandas>=2.0",
        "matplotlib>=3.7",
        "seaborn>=0.12",
        "biopython>=1.80",
        "scikit-learn>=1.2",
        "statsmodels>=0.14",
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov", "jupyter", "plotly", "kaleido"],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
