import numpy as np

from collections import defaultdict
from collections.abc import Iterable


class Energy():
    def __init__(self, data: np.array, info: dict = None, units: dict = None):
        if not isinstance(data, Iterable):
            raise TypeError("data has to be iterable.")

        if len(np.shape(data)) == 1:
            data = np.array([data])
        elif len(np.shape(data)) > 2:
            raise ValueError("data has to be a 1D or 2D array.")

        self.data = data

        self.__setup_info_dictionary__(info, units)

        # self.__make_attributes__()

    def __setup_info_dictionary__(self, info: dict = None, units: dict = None):
        """
        Sets up the info dictionary.

        If no info dictionary is given, a default dictionary is created, where
        the keys are the indices of the data array and the values are the indices
        of the data array. Furthermore a units dictionary can be given, where the
        where the keys have to match the keys of the info dictionary and the values
        are the units of the physical properties. If no units dictionary is given,
        the units are set to None.

        Parameters
        ----------
        info : dict
            A dictionary with the same length as the data array. The keys are either the
            indices of the data array or the names of the physical properties found in the 
            info file. The values are the indices of the data array.
        units : dict
            A dictionary with the same length as the data array. The keys are either the
            indices of the data array or the names of the physical properties found in the 
            info file. The values are the units of the physical properties.

        Raises
        ------
        TypeError
            If info is not a dictionary or None.
        TypeError
            if units is not a dictionary or None.
        ValueError
            If the length of info is not equal to the length of data.
        ValueError
            If the length of units is not equal to the length of data.
        ValueError
            If the keys of the info and units dictionary do not match.
        """
        if info is not None and not isinstance(info, dict):
            raise TypeError("info has to be a dictionary.")

        if units is not None and not isinstance(units, dict):
            raise TypeError("units has to be a dictionary.")

        if info is None:
            self.info_given = False
            info = defaultdict(lambda: None)
        else:
            self.info_given = True
            if len(info) != len(self.data):
                raise ValueError(
                    "The length of info dictionary has to be equal to the length of data.")

        if units is None:
            self.units_given = False
            units = defaultdict(lambda: None)
        else:
            self.units_given = True
            if len(units) != len(self.data):
                raise ValueError(
                    "The length of units dictionary has to be equal to the length of data.")

        self.info = info
        self.units = units

        if units.keys() != info.keys():
            raise ValueError(
                "The keys of the info and units dictionary have to be the same.")

        return info
