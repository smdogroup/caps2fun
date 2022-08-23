"""
Sean Engelstad
Georgia Tech SMDO
08/18/2022
"""

from typing import TYPE_CHECKING
from pathlib import Path

class FileManager:
    """
    Wrapper class for path and file management
    documentation of pathlib at https://docs.python.org/3/library/pathlib.html
    """
    def __init__(self):
        self._file_path = file

    @property
    def parent(self) -> Path:
        return self._file_path.root