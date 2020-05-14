from setuptools import setup


setup(
    name='superseal',
    python_requires='>3.6.6',
    version='0.0.1',
    url='https://github.com/stephenshank/superseal',
    download_url="https://github.com/stephenshank/superseal/archive/v0.0.1.tar.gz",
    description='Reference-guided viral quasipsecies reconstruction',
    author='Stephen D. Shank',
    author_email='sshank314@gmail.com',
    maintainer='Stephen D. Shank',
    maintainer_email='sshank314@gmail.com',
    install_requires=[
        'biopython>=1.73',
        'pandas>=0.24.2',
        'scipy>=1.1.0',
        'pysam>=0.15.2',
        'networkx>=2.3',
        'numpy>=1.16.4',
    ],
    packages=['superseal'],
    entry_points={
        'console_scripts': [
            'superseal = superseal.cli:command_line_interface'
        ]
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ],
    include_package_data=True
)
