def isinstance_of_list(obj, type):
    return isinstance(obj, list) and all(isinstance(o, type) for o in obj)
