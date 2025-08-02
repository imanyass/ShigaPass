#!/usr/bin/env python3

import os
import sys

from .blast import Blast
from .typing import Typing

class OptionsParser:
    """Class for options/arguments parser"""
    def __init__(self, version):
        self.version = version
        self._check_python()

    def _check_python(self):
        if sys.version_info.major < 3:
            print('Python 2 is no longer supported.')
            sys.exit(1)

    def parse_options(self, options):
        blast = Blast(options.db, options.threads)
        typing = Typing(blast)
        if not os.path.exists(options.outdir):
            os.makedirs(options.outdir)

        if options.mkdb:
            blast.make_blast_db()

        typing.run(options.list_file, options.outdir, options.keep_files)
