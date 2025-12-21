from amprenta_rag import constants


def test_constants_values():
    assert constants.DEFAULT_BATCH_SIZE == 100
    assert constants.DEFAULT_TOP_K == 10
    assert constants.MAX_CHUNKS == 50
    assert constants.DEFAULT_TTL == 3600
    assert constants.DEFAULT_PAGE_SIZE == 200

