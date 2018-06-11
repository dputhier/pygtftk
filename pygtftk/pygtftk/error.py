import pygtftk
from pygtftk.utils import message


class GTFtkError(Exception):

    def __init__(self, value):
        self.value = value

        if '__NON_INTERACTIVE__' in dir(pygtftk):

            message(value,
                    type="ERROR")
        else:
            raise GTFtkInteractiveError(value)

    def __str__(self):
        return repr(self.value)


class GTFtkInteractiveError(Exception):
    pass
