from __future__ import annotations

import logging
import textwrap
import os
import sys
import PQAnalysis.config as config

from pathlib import Path

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
                print_header(file=file)

    logger.propagate = False

    return logger


class CustomLogger(logging.Logger):
    def _log(self, level, msg, args, exc_info=None, extra=None, stack_info=False, stacklevel=1):
        if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
            stack_info = True

        super(CustomLogger, self)._log(
            level,
            msg,
            args,
            exc_info,
            extra,
            stack_info=stack_info,
            stacklevel=stacklevel
        )


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
