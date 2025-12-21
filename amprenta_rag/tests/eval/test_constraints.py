from amprenta_rag.eval import constraints


def test_no_pii_constraint_returns_empty():
    c = constraints.NoPIIConstraint()
    assert c.check("any text") == []


def test_role_header_constraint_returns_empty():
    c = constraints.RoleHeaderConstraint()
    assert c.check("ignored") == []


def test_constraint_violation_dataclass():
    violation = constraints.ConstraintViolation(code="X", message="m")
    assert violation.code == "X"
    assert violation.message == "m"

