from setuptools import setup
setup(name='coupledlogistic',
version='0.01',
description='Coupled Logistic Map',
url='https://github.com/artvalencio/coupled-logistic',
author='Arthur Valencio, IC-Unicamp, RIDC NeuroMat',
author_email='arthur@liv.ic.unicamp.br',
license='MIT',
packages=['coupledlogistic'],
classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
python_requires='>=3.6',
install_requires=['pandas','numpy'],
zip_safe=False)
