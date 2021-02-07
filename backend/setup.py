#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages


requirements = ['Bio','typing','pathlib',"graphviz","typing"]


test_requirements = ['pytest>=3']

setup(
    author="Group 5 (Dhruv Rathod, Ilya Yalchyk, Marlo Bareth, Pragya Mishra)",
    author_email='rathodhruv007@gmail.com, yalchik.ilya@gmail.com, marlo.bareth@gmail.com, pragyamishra004@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="The primary goal of this project is to create a pipeline that generates a 'co-author' network which "
                "details how strongly connected the queried author is to his/her co-authors.",
    entry_points={
        'console_scripts': [
            'authormaps=authormaps.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='authormaps',
    name='authormaps',
    packages=find_packages(include=['authormaps', 'authormaps.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
