import setuptools

setuptools.setup(
    name='alignpy',
    version=version,
    description='Align two FITS files to within a pixel',
    long_description=readme(),
    author='Annika Salmi',
    author_email='annika.salmi@yale.edu',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={"alignpy": ["data/*.pkl", "data/*.txt"]},
    install_requires=["numpy","scipy"],
    include_package_data=True,
    url='https://github.com/StarAlignment/alignpy',
    install_requires=[
        'numpy>=1.17',
        'scipy>=1',
        'astropy>=4'
        'python>3'
     ],
     classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
      ],
    python_requires='>=3.6',
)