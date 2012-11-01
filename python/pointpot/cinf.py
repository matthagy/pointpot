from hlab.pathutils import FilePath

import os
import ctypes
import warnings

from .paths import lib_dir

ctypes.pythonapi.PyFile_AsFile.argtypes = [ctypes.py_object]
ctypes.pythonapi.PyFile_AsFile.restype = ctypes.c_void_p


def select_lib(base_name, parent=None):
    if parent is None:
        parent = lib_dir
    base_path = parent.child(base_name)
    paths = list(path for path in (FilePath(base_path + suffix + ext)
                                   for suffix in ['', '_gcc']
                                   for ext in ['', '.dylib', '.so'])
                 if path.exists())
    if not paths:
        raise RuntimeError("couldn't locate path for base_name %r" % (base_name,))
    if len(paths) > 1:
        warnings.warn("located multiple libraries for base_name %r" % (base_name,), Warning)
    return paths[0]

def maybe_wire_dylib(basename, parent=None):
    if select_lib(basename, parent).endswith('.dylib'):
        wire_dylib()

wired_dylib = False
def wire_dylib():
    global wired_dylib
    if wired_dylib:
        return
    os.environ['DYLD_LIBRARY_PATH'] += ':' + lib_dir
    wired_dylib = True
