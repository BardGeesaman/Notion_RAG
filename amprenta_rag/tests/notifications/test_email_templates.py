from __future__ import annotations

import pytest
from amprenta_rag.notifications import email_templates as et

class FakeExperiment:
    def __init__(self, name="Exp", description="Desc", design_type="Design", organism="Human"):
        self.name = name
        self.description = description
        self.design_type = design_type
        self.organism = organism

class FakeDataset:
    def __init__(self, organism="Human"):
        self.organism = organism

def test_get_experiment_summary_html_basic():
    exp = FakeExperiment()
    html = et.get_experiment_summary_html(exp)
    assert "Experiment Summary" in html
    assert "Experiment Name:" in html
    assert "Exp" in html
    assert "Human" in html
    assert "0 dataset(s) associated" in html

def test_get_experiment_summary_html_with_datasets():
    exp = FakeExperiment(organism=None)
    ds1 = FakeDataset(organism="Mouse")
    ds2 = FakeDataset(organism="Rat")
    
    html = et.get_experiment_summary_html(exp, datasets=[ds1, ds2])
    assert "Mouse, Rat" in html
    assert "2 dataset(s) associated" in html

def test_get_experiment_summary_html_missing_fields():
    exp = FakeExperiment(name=None, description=None, design_type=None, organism=None)
    html = et.get_experiment_summary_html(exp)
    assert "Unknown Experiment" in html
    assert "No description available" in html
    assert "Unknown" in html

def test_get_digest_html():
    activities = [
        {"type": "experiment", "name": "Exp1", "url": "url1"},
        {"type": "experiment", "name": "Exp2", "url": "url2"},
        {"type": "compound", "name": "Cmp1", "url": "url3"},
    ]
    
    html = et.get_digest_html(activities, "Daily")
    assert "Daily Activity Digest" in html
    assert "Experiment (2)" in html
    assert "Compound (1)" in html
    assert "Exp1" in html
    assert "Cmp1" in html

def test_get_digest_html_truncated():
    activities = [{"type": "t", "name": f"Item {i}", "url": "#"} for i in range(15)]
    html = et.get_digest_html(activities, "Weekly")
    assert "... and 5 more" in html

def test_get_share_email_html():
    html = et.get_share_email_html("experiment", "Exp Name", "User A", "Check this out")
    assert "Shared Experiment: Exp Name" in html  # Header check might fail if text is different, let's check content
    assert "User A has shared a Experiment with you" in html
    assert "Exp Name" in html
    assert "Message:" in html
    assert "Check this out" in html

def test_get_share_email_html_no_message():
    html = et.get_share_email_html("compound", "Cmp Name", "User B")
    assert "Shared Compound" in html
    assert "Message:" not in html

