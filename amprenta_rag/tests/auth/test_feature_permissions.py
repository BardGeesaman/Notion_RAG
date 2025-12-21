def test_initialize_default_permissions_creates_defaults():
    session = FakeSession([]) # count = 0
    fp.initialize_default_permissions(session)
    
    assert len(session.added) > 0
    assert session.committed is True
    
    # Check Admin gets all (all set to visible)
    admin_perms = [p for p in session.added if p.role == "admin"]
    assert len(admin_perms) == len(fake_dash.ALL_PAGES)
    
    # Check Researcher gets all (some visible, some hidden)
    # The function adds permissions for ALL pages for the researcher role,
    # setting is_visible=False for admin pages.
    researcher_perms = [p for p in session.added if p.role == "researcher"]
    assert len(researcher_perms) == len(fake_dash.ALL_PAGES)
    
    # Verify specific visibility
    admin_page_perm = next(p for p in researcher_perms if p.page_name == "Admin")
    assert admin_page_perm.is_visible is False
    
    home_page_perm = next(p for p in researcher_perms if p.page_name == "Home")
    assert home_page_perm.is_visible is True
