from __future__ import annotations

import logging
import textwrap
import os
import sys
import types

from pathlib import Path

import PQAnalysis.config as config

from PQAnalysis.utils import print_header


def setup_logger(logger):
    stream_handler = logging.StreamHandler(stream=sys.stderr)
    stream_handler.setFormatter(CustomColorFormatter())
    logger.addHandler(stream_handler)

    if config.use_log_file:
        file_handler = logging.FileHandler(config.log_file_name)
        file_handler.setFormatter(CustomFormatter())
        logger.addHandler(file_handler)

        if os.stat(config.log_file_name).st_size == 0:
            with open(config.log_file_name, 'a') as file:
                print_header(stream=file)

    logger.propagate = False

    return logger


class CustomLogger(logging.Logger):
    def _log(self, level, msg, args, exception=None, exc_info=None, extra=None, stack_info=False, stacklevel=1):
        self._original_log(
            level,
            msg,
            args,
            exc_info,
            extra,
            stack_info,
            stacklevel
        )

        if level in [logging.CRITICAL, logging.ERROR]:
            if self.isEnabledFor(logging.DEBUG):
                try:
                    if exception is not None:
                        raise exception
                    else:
                        raise Exception
                except:
                    traceback = sys.exc_info()[2]
                    back_frame = traceback.tb_frame.f_back

                back_tb = types.TracebackType(
                    tb_next=None,
                    tb_frame=back_frame,
                    tb_lasti=back_frame.f_lasti,
                    tb_lineno=back_frame.f_lineno
                )

                if exception is not None:
                    raise Exception(msg).with_traceback(back_tb)
                else:
                    raise exception(msg).with_traceback(back_tb)

            else:
                exit(1)

    def _original_log(self, level, msg, args, exc_info=None, extra=None, stack_info=False, stacklevel=1):
        super()._log(level, msg, args, exc_info, extra, stack_info, stacklevel)

    def error(self, msg, exception=None, *args, **kwargs):
        if self.isEnabledFor(logging.ERROR):
            self._log(logging.ERROR, msg, args, exception, **kwargs)

    def original_error(self, msg, *args, **kwargs):
        if self.isEnabledFor(logging.ERROR):
            self._original_log(logging.ERROR, msg, args, **kwargs)

    def critical(self, msg, exception=None, *args, **kwargs):
        if self.isEnabledFor(logging.CRITICAL):
            self._log(logging.CRITICAL, msg, args, exception, **kwargs)

    def original_critical(self, msg, *args, **kwargs):
        if self.isEnabledFor(logging.CRITICAL):
            self._original_log(logging.CRITICAL, msg, args, **kwargs)


class CustomFormatter(logging.Formatter):

    level_keys = logging.getLevelNamesMapping().keys()
    longest_level = max(level_keys, key=len)

    def format_level(self, record):
        level = record.levelname + ":"
        additional_spaces = len(self.longest_level) - len(record.levelname)
        return level + ' ' * additional_spaces

    def format(self, record):
        message = record.msg.lstrip()
        record.msg = ''

        level = record.levelname
        formatted_level = self.format_level(record)

        level_keys = logging.getLevelNamesMapping().keys()
        longest_level_key = max(level_keys, key=len)

        name = record.name
        header = f'{formatted_level} {name}\n\n'

        messages = message.split('\n')
        wrapper = textwrap.TextWrapper(
            width=os.get_terminal_size().columns - len(level),
            initial_indent=' ' * (len(longest_level_key) + 2),
            subsequent_indent=' ' * (len(longest_level_key) + 2),
        )
        msg = '\n'.join(['\n'.join(wrapper.wrap(message))
                        for message in messages])

        record.msg = msg
        return header + msg + '\n'


class CustomColorFormatter(CustomFormatter):

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

    def format_level(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS["CUSTOM"])
        formatter = logging.Formatter(log_fmt)
        additional_spaces = len(self.longest_level) - len(record.levelname)
        return formatter.format(record) + ' ' * additional_spaces
