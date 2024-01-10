from beartype.typing import Callable
from decorator import decorator


@decorator
def extend_documentation(func, extension_function: Callable[[], str], *args, **kwargs):
    if extension_function.__doc__ is None:
        extension_function.__doc__ = ''

    if func.__doc__ is None:
        func.__doc__ = ''
    func.__doc__ += extension_function()

    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper
