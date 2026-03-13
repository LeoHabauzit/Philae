from setuptools import setup, find_packages

setup(
    name="philae",
    version="0.1.0",
    description="Description de votre projet",
    author="LeoHabauzit",
    url="https://github.com/LeoHabauzit/Philae",
    packages=find_packages(
        include=["Umat*", "simuEF*", "datas_simu*", "results_params*"]
    ),
    install_requires=[
        "numpy",
        "pandas",
        "IPython",
    ],
    dependency_links=[
        "git+https://github.com/simcoon/simcoon.git#egg=simcoon",
        "git+https://github.com/fedoo/fedoo.git#egg=fedoo",
    ],
    python_requires=">=3.8",
)
