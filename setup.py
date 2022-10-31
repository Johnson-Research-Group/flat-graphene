import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='flatgraphene',  
    version='0.3.6',
    scripts=['bin/flatgraphene'] ,
    author="Gabriel H. Brown",
    author_email="gabriel.h.brown@gmail.com",
    description="Generate mono and multilayer graphene geometries",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Johnson-Research-Group/flat-graphene",
    packages=setuptools.find_packages(),
    include_package_data = True, #include non-Python files specified in MANIFEST.in
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    keywords=[
        "twisted",
        "bilayer",
        "graphene",
        "graphite",
        "scientific",
        "engineering",
        "molecular dynamics",
        "atomistic",
    ],
    install_requires=[
        "numpy",
        "ase",
    ],
 )
