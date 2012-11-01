
import os
import sys
import distutils.sysconfig

import SCons
from SCons.Script.SConscript import SConsEnvironment


LIBPATH = ['/usr/lib', '.', '/usr/local/lib']
CPPPATH = ['/usr/local/include', '/usr/include']
CCFLAGS = ['-O3']

libname_base = 'pointpot'
# Some settings depending on the platform.
if sys.platform == 'darwin':
    libname = libname_base + '.dylib'
    linkflags = '-Wno-long-double -undefined suppress -flat_namespace '
    frameworksflags = '-flat_namespace -undefined suppress'
elif sys.platform == 'linux2':
    libname = libname_base + '.so'
    frameworksflags = ''
    linkflags = ''
else:
    raise SystemExit('Cannot build on %s.' % sys.platform)

sources = Glob('src/*.cpp')
includes = Glob('src/*.h')

AddOption('--prefix',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix')

env = Environment(LIBS=['m'], CPPPATH=CPPPATH, LIBPATH=LIBPATH,
                  SHLIBPREFIX='', CCFLAGS=CCFLAGS,
                  PREFIX = GetOption('prefix'))

SConsEnvironment.Chmod = SCons.Action.ActionFactory(os.chmod,
        lambda dest, mode: 'Chmod("%s", 0%o)' % (dest, mode))

SConsEnvironment.InstallHeader = lambda env, dest, files: InstallPerm(env, dest, files, 0644)

def InstallPerm(env, dest, files, perm):
    obj = env.Install(dest, files)
    for i in obj:
        env.AddPostAction(i, env.Chmod(str(i), perm))
    return dest

lib = env.SharedLibrary(libname, sources)

env.Alias('install-lib', env.Install('$PREFIX/lib', lib))
env.Alias('install-header', env.InstallHeader('$PREFIX/include', includes))

env.Alias('install', ['install-lib', 'install-header'])
