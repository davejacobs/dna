# Definitions for standard functions from functional programming

# 1. Currying

# Code to curry functions a la Ruby or Lisp
# Class taken from http://code.activestate.com/recipes/
# 52549-curry-associating-parameters-with-a-function/

class curry:
    def __init__(self, fun, *args, **kwargs):
        self.fun = fun
        self.pending = args[:]
        self.kwargs = kwargs.copy()

    def __call__(self, *args, **kwargs):
        a = self.pending + args
        self.kwargs.update(kwargs)
        return self.fun(*a, **self.kwargs)

# 2. Zipmap

def zipmap(f, coll):
    return dict(zip(coll, map(f, coll)))

# 3. Functional interpolation

def interpolate(string, info_map):
    return string % info_map
