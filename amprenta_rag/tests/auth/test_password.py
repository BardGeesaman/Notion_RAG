
from amprenta_rag.auth.password import hash_password, verify_password


def test_hash_password_produces_hash():
    plain = "correct horse battery staple"
    hashed = hash_password(plain)
    assert isinstance(hashed, str)
    assert hashed
    assert hashed != plain


def test_hash_password_different_each_time():
    plain = "same-password"
    h1 = hash_password(plain)
    h2 = hash_password(plain)
    assert h1 != h2


def test_verify_password_success():
    plain = "p@ssw0rd"
    hashed = hash_password(plain)
    assert verify_password(plain, hashed) is True


def test_verify_password_failure():
    plain = "p@ssw0rd"
    hashed = hash_password(plain)
    assert verify_password("wrong", hashed) is False


def test_verify_password_invalid_hash():
    # Should return False and not crash if the stored hash is malformed
    assert verify_password("anything", "not-a-valid-bcrypt-hash") is False


