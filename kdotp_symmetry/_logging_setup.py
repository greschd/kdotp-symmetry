# © 2017-2018, ETH Zurich, Institut für Theoretische Physik
# Author:  Dominik Gresch <greschd@gmx.ch>
"""
Creates the default kdotp-symmetry logger and logging handler.
"""

import sys
import logging

__all__ = ['LOGGER', 'DEFAULT_HANDLER']

LOGGER = logging.getLogger('kdotp_symmetry')
LOGGER.setLevel(logging.WARNING)

DEFAULT_HANDLER = logging.StreamHandler(sys.stdout)
LOGGER.addHandler(DEFAULT_HANDLER)
