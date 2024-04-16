"""
A module containing different decorators which could be useful.
"""

import time

from decorator import decorator
from collections import defaultdict
from inspect import getfullargspec
from beartype.typing import Any


def count_decorator(func, *args, **kwargs) -> Any:
    """
    Decorator which counts the number of times a function is called.
    """
    def wrapper(*args, **kwargs):
        """
        A wrapper function to count the number of times a function is called.
        """
        wrapper.counter += 1
        reset_counter = kwargs.get('reset_counter')

        if reset_counter:
            wrapper.counter = 1

        return func(*args, **kwargs)

    wrapper.counter = 0

    return wrapper


@decorator
def instance_function_count_decorator(func, *args, **kwargs):
    """
    Decorator which counts the number of times a function is called for an instance of a class.
    """
    self = args[0]

    if not hasattr(self, 'counter'):
        setattr(self, 'counter', defaultdict(int))
    self.counter[func.__name__] += 1

    reset_counter = get_argvar_by_name(func, args, 'reset_counter')

    if reset_counter:
        self.counter[func.__name__] = 1

    return func(*args, **kwargs)


@decorator
def timeit_in_class(func, *args, **kwargs) -> Any:
    """
    Decorator which measures the time a function of a class takes to execute
    and sets elapsed time as an attribute of the class.
    """
    def wrapper(*args, **kwargs) -> Any:
        """
        A wrapper function to measure the time a function takes to execute.

        Returns
        -------
        _type_
            _description_
        """
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        setattr(args[0], 'elapsed_time', end - start)
        return result

    return wrapper(*args, **kwargs)


def get_argvar_by_name(func, args, arg_name):
    """
    Returns the value of the argument with the given name.
    """
    arg_index = getfullargspec(func).args.index(arg_name)
    return args[arg_index]
