
__all__ = ["Block", "ForLoop"]

from ast import Param
from csmWriter.statements.base_statement import Statement
from csmWriter.parameters import *

class Block:
    def __init__(self, indent:int=0, main:bool=False, prefix=Statement(), suffix=Statement()):
        self._indent = indent
        self._is_main = main
        self._config_parameters = []
        self._design_parameters = []
        self._prefix = prefix
        self._statements = []
        self._suffix = suffix

    @property
    def indent(self) -> int:
        return self._indent

    @indent.setter
    def indent(self, new_indent:int):
        self._indent = new_indent
        for statement in self._statements:
            statement.indent = self._indent

    @property
    def statements(self):
        return self._statements

    @statements.setter
    def statements(self, new_statement):
        """
        add the new statement onto the list just by assignment process
        """
        if isinstance(new_statement, Parameter):
            assert(self._is_main)
            parameter = new_statement
            new_statement = parameter.statement
            if isinstance(parameter, ConfigParameter):
                self._config_parameters += [new_statement]
            elif isinstance(parameter, DesignParameter):
                self._design_parameters += [new_statement]
        elif isinstance(new_statement, Block):
            new_statement.indent = self.indent + 1
            self._statements += [new_statement]
        else:
            self._statements += [new_statement]

    def __str__(self):
        print(self._design_parameters)
        return "".join([str(param) for param in self._config_parameters]) + \
            "".join([str(param) for param in self._design_parameters]) + \
            str(self._prefix) + \
            "".join([str(statement) for statement in self._statements]) + \
            str(self._suffix)

class ForLoop(Block):
    def __init__(self, index_name:str="foo", num:int=1):
        prefix=Statement(content=f"patbeg  {index_name}  {num}")
        suffix=Statement(content="patend")
        super(ForLoop,self).__init__(prefix=prefix, suffix=suffix)
        self._index = Parameter(name=index_name)
        
    @property
    def index(self) -> Parameter:
        return self._index
