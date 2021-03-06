import setuptools

setuptools.setup(
    name="igemtools_lib",
    version="0.0.1",
    author="Matthijs Tadema",
    author_email="M.J.Tadema@gmail.com",
    description="Bioinformatics toolkit for the iGEM Groningen 2018 team",
    #long_description=long_description,
    long_description_content_type="text/markdown",
    #url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    scripts=[
        "bin/igemtools.py"
    ],
    install_requires=[
        'biopython',
        'numpy'
    ]
)
