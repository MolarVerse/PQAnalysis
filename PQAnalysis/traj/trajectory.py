import PQAnalysis.traj.frame as frame


def read(filename):
    '''
    Read a trajectory from a file.
    '''

    if filename.endswith('.xyz'):
        return read_xyz(filename)


def read_xyz(filename):
    '''
    Read a trajectory from a .xyz file.
    '''

    frames = []
    with open(filename, 'r') as f:
        frame_string = ''
        for line in f:
            if line.strip() == '':
                frame_string += line
                continue
            elif line.split()[0].isdigit():
                if frame_string != '':
                    frames.append(frame.read_frame(frame_string))
                frame_string = line
            else:
                frame_string += line
        frames.append(frame.read_frame(frame_string))

    return Trajectory(frames)


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
