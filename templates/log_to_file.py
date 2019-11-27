#!/usr/bin/env python
"""Logging setup suitable for cronjobs.

A file gets the detailed logs while stderr gets warnings, errors and critical.
"""

import logging
import os
import sys

def main():
    logger.info("First entry as INFO")
    logger.warning("Second entry as WARNING")

if __name__ == '__main__':
    logfile = "logging_test.txt"
    logfile_directory = os.path.dirname(os.path.abspath(logfile))
    if not os.path.exists(logfile_directory):
        logging.error("The directory for the specified log file does not exist. Aborting")
        sys.exit(-1)

    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # Handler that will log warnings or worse
    stderr_handler = logging.StreamHandler()
    stderr_handler.setLevel(logging.WARNING)
    logger.addHandler(stderr_handler)

    main()
