"""
A module containing different decorators which could be useful.
"""


def count_decorator(func):
    """
    Decorator which counts the number of times a function is called.
    """
    def wrapper(*args, **kwargs):
        wrapper.counter += 1
        return func(*args, **kwargs)
    wrapper.counter = 0
    return wrapper
