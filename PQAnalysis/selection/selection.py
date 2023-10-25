from PQAnalysis.utils.instance import isinstance_of_list

class Selection:
    def __init__(self, selection):
        self.selection = selection


class AtomSelection(Selection):
    def __init__(self, selection, frame=None):
        if isinstance(selection, int):
            super().__init__([selection])
        elif isinstance_of_list(selection, int):
            super().__init__(selection)
        elif frame is None:
            raise ValueError('Frame must be provided when selection is not an integer.')
        elif isinstance(selection, str):


        
    #     if 
    #     elif isinstance(selection, str):
    #         if frame is None:
    #             raise ValueError('Frame must be provided when selection is a string.')
    #         else:

    # def parse_selection_string(self, selection_string, frame):
    #     frame.atoms.index(selection_string)


