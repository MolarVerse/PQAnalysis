import logging
import textwrap
import os


class CustomFormatter(logging.Formatter):

    bold_yellow = "\x1b[33;1m"
    bold_green = "\x1b[32;1m"
    bold_orange = "\x1b[33;1m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    header = "%(levelname)s -"

    FORMATS = {
        logging.DEBUG: bold_yellow + header + reset,
        logging.INFO: bold_green + header + reset,
        logging.WARNING: bold_orange + header + reset,
        logging.ERROR: red + header + reset,
        logging.CRITICAL: bold_red + header + reset
    }

    """Multi-line formatter."""

    def _format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

    def format(self, record):
        message = record.msg.lstrip()
        record.msg = ''
        level = record.levelname
        color_level = self._format(record)
        name = record.name
        header = f'{color_level} {name}:\n\n'
        wrapper = textwrap.TextWrapper(
            width=os.get_terminal_size().columns - 10,
            initial_indent='+++' + ' ' * len(level),
            subsequent_indent='+++' + ' ' * len(level)
        )
        msg = wrapper.wrap(message)
        record.msg = msg
        return header + "\n".join(msg) + '\n'
