

class AutoRepr(object):

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,
                           ', '.join(map(repr, self.repr_args())))

    def __str__(self):
        return repr(self)


