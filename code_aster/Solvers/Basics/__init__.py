# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

"""
Definition of common objects for non linear operators.
"""

__all__ = (
    "DispatcherMixin",
    "EventId",
    "EventSource",
    "Observer",
    "ProblemType",
    "Context",
    "ContextMixin",
    "LoggingManager",
    "PhysicalState",
    "Residuals",
)

from .bases import DispatcherMixin, EventId, EventSource, Observer, ProblemType
from .context import Context, ContextMixin
from .logging_manager import LoggingManager
from .physical_state import PhysicalState
from .residual import Residuals
