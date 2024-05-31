"""
A module containing custom logging classes and functions,
that are used in the PQAnalysis package.
"""

import logging

import textwrap
import os
import shutil
import sys
import types

from beartype.typing import Any

from PQAnalysis.config import log_file_name, use_log_file
from PQAnalysis.utils import print_header
from PQAnalysis.exceptions import PQException



def setup_logger(logger: logging.Logger) -> logging.Logger:
    """
    general setup for the logger

    This function can be used to setup a generic logger for functions/classes...
    of the PQAnalysis package. It automatically adds a stream handler to the logger
    and if the use_log_file is set to True, it also adds a file handler to
    the logger. The file handler writes the log messages to the file specified in
    the config log_file_name.

    Parameters
    ----------
    logger : logging.Logger
        The logger to setup.

    Returns
    -------
    logging.Logger
        The logger that was setup.
    """
    if logging.getLogger(logger.name).handlers == []:
        stream_handler = logging.StreamHandler(stream=sys.__stderr__)
        stream_handler.setFormatter(CustomColorFormatter())
        logger.addHandler(stream_handler)

        if use_log_file:
            file_handler = logging.FileHandler(log_file_name)
            file_handler.setFormatter(CustomFormatter())
            logger.addHandler(file_handler)

            if os.stat(log_file_name).st_size == 0:
                with open(log_file_name, 'a', encoding='utf-8') as file:
                    print_header(file=file)
    else:
        logger = logging.getLogger(logger.name)

    logger.propagate = False

    return logger



class CustomLogger(logging.Logger):

    """
    A custom logger class that extends the logging.Logger class.

    This class extends the logging.Logger class and is set as the
    default logger class for the PQAnalysis package. It re-implements
    the _log method of the logging.Logger class in conjunction with
    the error and critical methods. For more details see the
    documentation of the _log, error and critical methods.

    To return to the original logging.Logger class, the original_error
    and original_critical methods can be used.
    """

    def _log(
        self,  # pylint: disable=arguments-differ
        level: Any,
        msg: Any,
        args: Any,
        exception: Exception | None = None,
        **kwargs
    ) -> None:
        """
        This method is a wrapper method for the original _log method of the logging.Logger class.

        It logs the message with the given level and raises an exception if the level 
        is logging.ERROR or logging.CRITICAL and the logger is enabled for logging.DEBUG.
        If the logger is not enabled for logging.DEBUG, the program exits with code 1.


        Parameters
        ----------
        level : Any
            The level of the log message.
        msg : Any
            The message to log.
        args : Any
            The arguments of the log message.
        exception : Exception, optional
            The exception to raise, by default None.

        Raises
        ------
        exception
            if the level is logging.ERROR or logging.CRITICAL and the logger
            is enabled for logging.DEBUG.
        Exception
            if the level is logging.ERROR or logging.CRITICAL and the logger
            is not enabled for logging.DEBUG.
        """
        if exception is not None:
            extra = {'custom_exception': exception}
        else:
            extra = None

        if 'extra' in kwargs:
            del kwargs['extra']

        self._original_log(level, msg, args, extra, **kwargs)

        if level in [logging.CRITICAL, logging.ERROR]:

            exception = exception or Exception

            if self.isEnabledFor(logging.DEBUG):
                back_tb = None

                try:
                    raise exception  # pylint: disable=broad-exception-raised
                except exception:  # pylint: disable=broad-except
                    traceback = sys.exc_info()[2]
                    back_frame = traceback.tb_frame.f_back

                    back_tb = types.TracebackType(
                        tb_next=None,
                        tb_frame=back_frame,
                        tb_lasti=back_frame.f_lasti,
                        tb_lineno=back_frame.f_lineno
                    )

                raise exception(msg).with_traceback(back_tb)

            def exception_hook(exc_type, exc_value, exc_traceback):
                """
                A custom exception hook that ignores the CustomLoggerException.
                """
                if isinstance(exc_type, PQException):
                    sys.__excepthook__(exc_type, exc_value, exc_traceback)

            sys.excepthook = exception_hook
            raise exception(msg)  # pylint: disable=broad-exception-raised

    def _original_log(
        self, level: Any, msg: Any, args: Any, extra=None, **kwargs
    ) -> None:
        """
        The original _log method of the logging.Logger class.

        This method logs the message with the given level and
        the given arguments and keyword arguments.

        Parameters
        ----------
        level : Any
            The level of the log message.
        msg : Any
            The message to log.
        args : Any
            The arguments of the log message.
        extra : Any, optional
            Additional information to log, by default None.
        **kwargs : Any
            The keyword arguments of the log message.
        """

        super()._log(level, msg, args, extra=extra, **kwargs)

    def error(
        self,
        msg: Any,
        *args,
        exception: Exception | None = None,
        **kwargs
    ) -> None:
        """
        This method logs the message with the logging.ERROR level and raises
        an exception if the logger is enabled for logging.DEBUG. If the logger
        is not enabled for logging.DEBUG, the program exits with code 1. If 
        an exception is given, it is raised with the message of the log message
        otherwise a generic Exception is raised.

        In order to return to the original logging.Logger class, the original_error
        method can be used.

        Parameters
        ----------
        msg : Any
            The message to log.
        exception : Exception, optional
            The exception to raise, by default None. If None, a generic Exception is raised.
        *args : Any
            The arguments of the log message.
        **kwargs : Any
            The keyword arguments of the log message.
        """
        if self.isEnabledFor(logging.ERROR):
            self._log(logging.ERROR, msg, args, exception=exception, **kwargs)

    def original_error(self, msg: Any, *args, **kwargs) -> None:
        """
        This method logs the message with the logging.ERROR level.

        It corresponds to the original error method of the logging.Logger class.

        Parameters
        ----------
        msg : Any
            The message to log.
        *args : Any
            The arguments of the log message.
        **kwargs : Any
            The keyword arguments of the log message.
        """

        if self.isEnabledFor(logging.ERROR):
            self._original_log(logging.ERROR, msg, args, **kwargs)

    def critical(
        self,
        msg: Any,
        *args,
        exception: Exception | None = None,
        **kwargs
    ) -> None:
        """
        This method logs the message with the logging.CRITICAL level and
        raises an exception if the logger is enabled for logging.DEBUG. 
        If the logger is not enabled for logging.DEBUG, the program exits
        with code 1. If an exception is given, it is raised with the 
        message of the log message otherwise a generic Exception is raised.

        Parameters
        ----------
        msg : Any
            The message to log.
        exception : Exception, optional
            The exception to raise, by default None. If None, a generic Exception is raised.
        *args : Any
            The arguments of the log message.
        **kwargs : Any
            The keyword arguments of the log message.
        """
        if self.isEnabledFor(logging.CRITICAL):
            self._log(
                logging.CRITICAL,
                msg,
                args=args,
                exception=exception,
                **kwargs
            )

    def original_critical(self, msg: Any, *args, **kwargs) -> None:
        """
        This method logs the message with the logging.CRITICAL level.

        It corresponds to the original critical method of the logging.Logger class.

        Parameters
        ----------
        msg : Any
            The message to log.
        *args : Any
            The arguments of the log message.
        **kwargs : Any
            The keyword arguments of the log message.
        """
        if self.isEnabledFor(logging.CRITICAL):
            self._original_log(logging.CRITICAL, msg, args, **kwargs)



class CustomFormatter(logging.Formatter):

    """
    A custom formatter class that extends the logging.Formatter class.

    This class extends the logging.Formatter class and is used to format
    the log messages of the PQAnalysis package. It re-implements the 
    format method of the logging.Formatter class to format the log 
    messages in a custom way. The format method is used to format the log
    messages of the logger.
    """
    level_keys = logging.getLevelNamesMapping().keys()
    longest_level = max(level_keys, key=len)

    def format_level(self, record: logging.LogRecord) -> str:
        """
        Formats the level of the log message.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to format.

        Returns
        -------
        str
            The formatted level of the log message.
        """
        level = record.levelname + ":"
        additional_spaces = len(self.longest_level) - len(record.levelname)
        return level + ' ' * additional_spaces

    def format(self, record: logging.LogRecord) -> str:
        """
        Formats the log message.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to format.

        Returns
        -------
        str
            The formatted log message.
        """

        message = record.msg.lstrip()
        record.msg = ''

        level = record.levelname
        formatted_level = self.format_level(record)

        level_keys = logging.getLevelNamesMapping().keys()
        longest_level_key = max(level_keys, key=len)

        name = record.name

        if hasattr(record, 'custom_exception'):
            header = (
                f'{formatted_level} {name} - '
                f'{record.custom_exception.__qualname__}\n\n'
            )
        else:
            header = f'{formatted_level} {name}\n\n'

        messages = message.split('\n')
        wrapper = textwrap.TextWrapper(
            width=shutil.get_terminal_size(fallback=(80, 100)).columns -
            len(level),
            initial_indent=' ' * (len(longest_level_key) + 2),
            subsequent_indent=' ' * (len(longest_level_key) + 2),
        )
        msg = '\n'.join(
            ['\n'.join(wrapper.wrap(message)) for message in messages]
        )

        record.msg = msg
        return '\n' + header + msg



class CustomColorFormatter(CustomFormatter):

    """
    A custom color formatter class that extends the CustomFormatter class.

    This class extends the CustomFormatter class, so that for each log level
    a different color is used. The colors are defined in the FORMATS dictionary.
    The color is set to the log level of the log message.
    """

    bold_yellow = "\x1b[33;1m"
    bold_green = "\x1b[32;1m"
    bold_orange = "\x1b[33;1m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    bold_cyan = "\x1b[36;1m"
    header = "%(levelname)s:"

    FORMATS = {
        logging.DEBUG: bold_yellow + header + reset,
        logging.INFO: bold_green + header + reset,
        logging.WARNING: bold_orange + header + reset,
        logging.ERROR: red + header + reset,
        logging.CRITICAL: bold_red + header + reset,
        "CUSTOM": bold_cyan + header + reset
    }

    def format_level(self, record: logging.LogRecord) -> str:
        """
        Formats the level of the log message.

        It formats the level of the log message with the color of the log level.

        Parameters
        ----------
        record : logging.LogRecord
            the log record to format

        Returns
        -------
        str
            the formatted level of the log message
        """

        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS["CUSTOM"])
        formatter = logging.Formatter(log_fmt)
        additional_spaces = len(self.longest_level) - len(record.levelname)
        return formatter.format(record) + ' ' * additional_spaces
