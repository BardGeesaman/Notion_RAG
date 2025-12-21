from __future__ import annotations

from uuid import uuid4
import importlib
import sys
from types import ModuleType

import pytest


def _install_fake_pptx():
    fake_pptx = ModuleType("pptx")
    fake_pptx.util = ModuleType("pptx.util")
    fake_pptx.util.Pt = lambda x: x
    fake_pptx.Presentation = None  # will be set later
    sys.modules["pptx"] = fake_pptx
    sys.modules["pptx.util"] = fake_pptx.util
    return fake_pptx


class _FakeTextFrame:
    def __init__(self):
        self.paragraphs = []

    def clear(self):
        self.paragraphs.clear()

    def add_paragraph(self):
        p = type("P", (), {"text": "", "level": 0, "font": type("F", (), {"size": None})()})()
        self.paragraphs.append(p)
        return p


class _FakePlaceholder:
    def __init__(self):
        self.text_frame = _FakeTextFrame()
        self.text = ""


class _FakeSlide:
    def __init__(self):
        self.shapes = type("S", (), {"title": type("T", (), {"text": ""})()})()
        self.placeholders = [None, _FakePlaceholder()]


class _FakeSlides:
    def __init__(self):
        self.slides = []

    def add_slide(self, layout):
        slide = _FakeSlide()
        self.slides.append(slide)
        return slide


class _FakePresentation:
    def __init__(self):
        self.slide_layouts = [0, 1]
        self.slides = _FakeSlides()

    def save(self, output):
        output.write(b"fake")


class _FakeQuery:
    def __init__(self, obj):
        self.obj = obj

    def filter(self, *_):
        return self

    def first(self):
        return self.obj


class _FakeDB:
    def __init__(self, obj):
        self.obj = obj

    def query(self, model):
        return _FakeQuery(self.obj)


def test_generate_experiment_slides(monkeypatch):
    fake_pptx = _install_fake_pptx()
    import amprenta_rag.export.slide_generator as sg

    monkeypatch.setattr(sg, "Presentation", _FakePresentation)
    monkeypatch.setattr(sg, "Pt", lambda x: x)
    fake_pptx.Presentation = _FakePresentation

    class FakeDataset:
        def __init__(self, organism=None):
            self.organism = organism

    class FakeExperiment:
        def __init__(self):
            self.id = uuid4()
            self.name = "Exp"
            self.created_at = None
            self.description = "desc"
            self.design_type = "design"
            self.datasets = [FakeDataset("human")]
            self.disease = ["dis"]
            self.matrix = ["m"]
            self.sample_groups = []

    db = _FakeDB(FakeExperiment())
    result = sg.generate_experiment_slides(uuid4(), db)
    assert isinstance(result, (bytes, bytearray))


def test_generate_dataset_slides(monkeypatch):
    fake_pptx = _install_fake_pptx()
    import amprenta_rag.export.slide_generator as sg

    monkeypatch.setattr(sg, "Presentation", _FakePresentation)
    monkeypatch.setattr(sg, "Pt", lambda x: x)
    fake_pptx.Presentation = _FakePresentation

    class FakeDataset:
        def __init__(self):
            self.id = uuid4()
            self.name = "DS"
            self.created_at = None
            self.omics_type = "omics"
            self.dataset_source_type = "src"
            self.description = "desc"
            self.organism = ["human"]
            self.disease = ["d"]
            self.sample_type = ["blood"]
            self.data_origin = None
            self.summary = None
            self.sample_groups = []

    db = _FakeDB(FakeDataset())
    result = sg.generate_dataset_slides(uuid4(), db)
    assert isinstance(result, (bytes, bytearray))


def test_generate_experiment_slides_missing(monkeypatch):
    fake_pptx = _install_fake_pptx()
    import amprenta_rag.export.slide_generator as sg

    monkeypatch.setattr(sg, "Presentation", _FakePresentation)
    monkeypatch.setattr(sg, "Pt", lambda x: x)
    fake_pptx.Presentation = _FakePresentation

    db = _FakeDB(None)
    with pytest.raises(ValueError):
        sg.generate_experiment_slides(uuid4(), db)

