import sys
import os


class BaseWriter:

    def __init__(self, filename, mode):
        if mode not in ['w', 'a']:
            raise ValueError('Invalid mode has to be either \'w\' or \'a\'.')
        elif mode == 'w' and filename is not None and os.path.isfile(filename):
            raise ValueError(
                f"File {filename} already exists. Use mode \'a\' to append to file.")

        if filename is None:
            self.file = sys.stdout
        else:
            self.file = None

        self.mode = 'a'
        self.filename = filename
        self.format = format

    def open(self):
        if self.filename is not None:
            self.file = open(self.filename, self.mode)

    def close(self):
        if self.file is not None:
            self.file.close()
