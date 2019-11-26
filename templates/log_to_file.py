#!/usr/bin/env python
"""Logging setup suitable for cronjobs.

A file gets the detailed logs while stderr gets warnings, errors and critical.
"""

import logging

def main():
    logger.info("First entry as INFO")
    logger.warning("Second entry as WARNING")

if __name__ == '__main__':
    logging.basicConfig(filename="logging_test.txt", level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # Handler that will log warnings or worse
    stderr_handler = logging.StreamHandler()
    stderr_handler.setLevel(logging.WARNING)
    logger.addHandler(stderr_handler)

    main()
