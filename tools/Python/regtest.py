#!/usr/bin/env python3

"""
A simple regression test framework for a AMReX-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.
They are separated as such in this file.

This test framework understands source based out of the AMReX framework.

"""

from __future__ import print_function

from pathlib import Path
import asyncio
import sys

from rtl2 import regtest

if __name__ == "__main__":
    todo = regtest.main(sys.argv[1:])

    loop = asyncio.get_event_loop()
    n = loop.run_until_complete(todo)

    sys.exit(n)
