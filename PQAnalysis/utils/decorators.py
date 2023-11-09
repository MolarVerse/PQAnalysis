"""
A module containing different decorators which could be useful.
"""

from collections import defaultdict


def count_decorator(func):
    """
    Decorator which counts the number of times a function is called.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Returns
    -------
    wrapper : function
        Decorated function.
    """
    def wrapper(*args, **kwargs):
        wrapper.counter += 1
        reset_counter = kwargs.get('reset_counter')

        if reset_counter:
            wrapper.counter = 1

        return func(*args, **kwargs)
    wrapper.counter = 0
    return wrapper


def instance_function_count_decorator(func):
    """
    Decorator which counts the number of times a function is called for an instance of a class.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Returns
    -------
    new_func : function
        Decorated function.
    """

    def new_func(self, *args, **kwargs):
        if not hasattr(self, 'counter'):
            setattr(self, 'counter', defaultdict(int))
        self.counter[func.__name__] += 1

        reset_counter = kwargs.get('reset_counter')

        if reset_counter:
            self.counter[func.__name__] = 1

        func(self, *args, **kwargs)
    new_func.__name__ = func.__name__
    return new_func
