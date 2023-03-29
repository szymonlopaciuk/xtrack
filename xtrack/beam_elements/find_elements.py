# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #
import inspect

from types import ModuleType
from typing import Set

from xtrack import BeamElement


def find_elements_in_module(
        module: ModuleType,
        _visited: Set[ModuleType] = None,
) -> Set[BeamElement]:
    """Get all beam elements in a module and its submodules. Runs a simple
    breadth-first search.

    Args:
        module: The module to search.
        _visited: A set of modules to exclude from the search.

    Returns:
        A set of beam elements.
    """
    beam_elements = set()
    _visited = _visited or set()

    if module in _visited:
        return set()
    _visited.add(module)

    for name, obj in inspect.getmembers(module):
        if isinstance(obj, type) and issubclass(obj, BeamElement):
            beam_elements.add(obj)
        elif isinstance(obj, ModuleType):
            beam_elements.update(find_elements_in_module(obj, _visited))

    return beam_elements
