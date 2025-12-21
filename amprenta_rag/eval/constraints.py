from __future__ import annotations

from dataclasses import dataclass
from typing import List, Protocol


@dataclass(frozen=True)
class ConstraintViolation:
    code: str
    message: str


class OutputConstraint(Protocol):
    def check(self, output_text: str) -> List[ConstraintViolation]:
        ...


class NoPIIConstraint:
    """
    Very lightweight placeholder. Extend with proper checks.
    """

    def check(self, output_text: str) -> List[ConstraintViolation]:
        # Placeholder: no-op
        return []


class RoleHeaderConstraint:
    """
    Example format guard: if using role headers, ensure they are present.
    Disabled by default (opt-in in evaluator).
    """

    def check(self, output_text: str) -> List[ConstraintViolation]:
        del output_text  # unused in placeholder
        return []


