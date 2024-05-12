"""
A module containing different decorators which could be useful.
"""

import time

from collections import defaultdict
from inspect import getfullargspec

from decorator import decorator
from beartype.typing import Any



def count_decorator(func, *args, **kwargs) -> Any:  # pylint: disable=unused-argument, missing-function-docstring
    """
    Decorator which counts the number of times a function is called.

    Parameters
    ----------
    func : Callable
        The function to be decorated.
    args : 
        The arguments of the function.
    kwargs : 
        The keyword arguments of the function.
    """

    def wrapper(*args, **kwargs):
        """
        A wrapper function to count the number of times a function is called.
        """
        wrapper.counter += 1

        if kwargs.get('reset_counter', False):
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

    if get_arg_var_by_name(func, args, 'reset_counter'):
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
        """
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        setattr(args[0], 'elapsed_time', end - start)
        return result

    return wrapper(*args, **kwargs)



def get_arg_var_by_name(func, args, arg_name):
    """
    Returns the value of the argument with the given name.
    """
    arg_index = getfullargspec(func).args.index(arg_name)
    return args[arg_index]
