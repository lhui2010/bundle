"""
Function of this scripts
"""

import logging
import coloredlogs

from iga.apps.base import emain

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


def abc(foo=None, bar=None):
    """
    Remember:
     1. Positional args are recognized by software with None value
     2. All Non-positional args need to be given a default value, like ''
     2. no % is allowed in function docs.
    Returns:
    """

if __name__ == "__main__":
    emain()
