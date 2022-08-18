
__all__ = ["CsmFile"]

from hashlib import new
from typing import TYPE_CHECKING
from .block.base_block import Block

class CsmFile:
    def __init__(self, filename:str):
        assert(".csm" in filename)
        self._filename = filename
        self._main_block = Block(main=True)

    @property
    def main(self) -> Block:
        return self._main_block.statements

    @main.setter
    def main(self, new_statement):
        self._main_block.statements = new_statement

    def write(self):
        handle = open(self._filename, "w")
        handle.write(str(self._main_block))