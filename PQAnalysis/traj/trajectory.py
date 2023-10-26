import PQAnalysis.traj.frame as frame
import numpy as np


class Trajectory:
    '''
    A trajectory is a sequence of frames.
    '''

    def __init__(self, frames=None):
        if frames is None:
            self.frames = []
        else:
            self.frames = frames

    def __len__(self):
        return len(self.frames)

    def __getitem__(self, key):
        return self.frames[key]

    def __iter__(self):
        return iter(self.frames)

    def __contains__(self, item):
        return item in self.frames

    def __next__(self):
        return next(self.frames)

    def __add__(self, other):
        return Trajectory(self.frames + other.frames)
