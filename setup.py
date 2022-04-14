from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    author="Jiaying Lai",
    description="A package for predicting mutation pathogenicity",
    name="lyrus",
    version="0.1.0",
    packages=find_packages(include=["LYRUS","LYRUS.*"]),
    package_data={'LYRUS': ['train.csv']},
    install_requires=[
        'bio>=0.4.1',
        'numpy>=1.18.5',
        'beautifulsoup4>=4.7.1',
        'variation-number>=0.1.7'
    ],
    python_requires='>=3.7',
    long_description=long_description,
    long_description_content_type='text/markdown',
    include_package_data=True
)