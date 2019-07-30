from setuptools import setup

dev_dependencies = [
    'xlrd==1.2.0'
]

setup(
    name='convex_qsr',
    version='0.1.2',
    url='https://github.com/stephenshank/convex-qsr',
    download_url="https://github.com/stephenshank/convex-qsr/archive/v0.1.1.tar.gz",
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
        'flask>=1.1.1'
    ],
    tests_require=dev_dependencies,
    extras_require={
        'dev': dev_dependencies
    },
    packages=['convex_qsr'],
    entry_points={
        'console_scripts': [
            'cqsr_ec = convex_qsr.cli:error_correction',
            'cqsr_rg = convex_qsr.cli:read_graph',
            'cqsr_qr = convex_qsr.cli:quasispecies_reconstruction',
        ]
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ],
)
