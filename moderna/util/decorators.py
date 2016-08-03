#!/usr/bin/env python
"""
Decorators for moderna.py.
"""


from .Errors import ModernaError
from .LogFile import log
from moderna.examples.usage_examples import COMMAND_EXAMPLES


def simple_decorator(decorator):
    """This decorator can be used to turn simple functions
    into well-behaved decorators, so long as the decorators
    are fairly simple. If a decorator expects a function and
    returns a function (no descriptors), and if it doesn't
    modify function attributes or docstring, then it is
    eligible to use this. Simply apply @simple_decorator to
    your decorator and it will automatically preserve the
    docstring and function attributes of functions to which
    it is applied."""
    def new_decorator(f):
        g = decorator(f)
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        g.__dict__.update(f.__dict__)
        return g
    # Now a few lines needed to make simple_decorator itself
    # be a well-behaved decorator.
    new_decorator.__name__ = decorator.__name__
    new_decorator.__doc__ = decorator.__doc__
    new_decorator.__dict__.update(decorator.__dict__)
    return new_decorator

#
# decorator that catches ModernaErrors and writes them to log file.
#
@simple_decorator
def toplevel_function(func):
    """If a ModeRNA Error occurs, it will be written to the logfile.
    Also adds examples to the documentation strings.
    """
    if not func.__name__ in COMMAND_EXAMPLES:
        raise ValueError("No example for ModeRNA function: %s"%func.__name__)
    func.__doc__ += '\nExamples:\n%s' % COMMAND_EXAMPLES[func.__name__]

    def catch_and_log_errors(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ModernaError:
            log.write_error()
            if log.raise_exceptions:
                raise

    return catch_and_log_errors
