import setuptools
setuptools.setup(
    name="holteandtalley",
    version="0.0.2",
    author="Garrett Finucane",
    author_email="garrettdreyfus@gmail.com",
    description="A python adaptation of the Holte and Talley mixed layer depth algorithm",
    url="https://github.com/garrettdreyfus/python-holteandtalley",
    license = "MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy',
    ],
    packages=['holteandtalley'],
)
