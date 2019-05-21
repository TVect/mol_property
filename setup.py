# from distutils.core import setup
from setuptools import setup

setup(name="mol_property",
      version="1.0.0",
      description="Prediction of pKa from chemical structure using machine learning approaches",
      long_description=open('README.md', encoding='utf-8').read(),
      author="Chin",
      author_email="chin340823@163.com",
      url="https://github.com/TVect/mol_property",
      packages=['mol_property', 'mol_property.pka', 'mol_property.solubility'],
      package_data={'mol_property': ["mol_property/pka/model/*.pkl"]},
      include_package_data=True,
      #install_requires=['numpy>=1.16.3',
      #                  'pandas>=0.24.2',
      #                  'rdkit>=2009.Q1-1',
      #                  'scikit_learn==0.21.1'],
      classifiers=['Development Status :: 2 - Pre-Alpha', 
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3.6',
                   'Topic :: Scientific/Engineering :: Chemistry'],
      zip_safe=False,
) 

