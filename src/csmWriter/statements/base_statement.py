
__all__ = ["Statement"]

from typing import TYPE_CHECKING

class Statement:
    def __init__(self, content:str='', indent:int=0):
        """
        Statement base class, creates blank line
        """
        self._indent = indent
        self._tab = '\t'
        self._content = content

    @property
    def indent(self) -> int:
        return self._indent

    @indent.setter
    def indent(self, value):
        self._indent = value

    @property
    def content(self) -> str:
        return self._content

    def __str__(self):
        return self._tab * self._indent + self.content + '\n'