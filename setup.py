from setuptools import setup

setup(
    name="M2R-Gillespie-Algorithm",
    version="0.1.0",
    description="Gillespie SSA implementations for stochastic biochemical systems",
    author="Alex Jin, Yanbo Huang, Jiali Xu, Yuhan Zhu, Xuancheng Yi",
    license="MIT",
    package_dir={"": "src"},
    py_modules=["gillespie"],
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
    ],
    python_requires=">=3.7",
)
