"""
A package containing classes and functions to handle radial distribution functions.

The rdf package contains the following submodules:
    
        - rdf
        
The rdf package contains the following classes:

        - RadialDistributionFunction
        
The rdf package contains the following exceptions:

        - RDFError
        - RDFWarning
"""

from .exceptions import RDFError, RDFWarning
from .rdf import RadialDistributionFunction
