"""Tests for retrosynthesis service."""

import pytest
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.chemistry.retrosynthesis import (
    RetrosynthesisPlanner,
    score_route,
    check_building_blocks,
    SynthesisRoute,
    SynthesisStep,
    SynthesisTree,
    RouteScore,
)


class TestRetrosynthesisPlanner:
    def test_mock_planner_returns_valid_tree(self):
        """Mock planner returns SynthesisTree with routes."""
        planner = RetrosynthesisPlanner(backend="mock")
        tree = planner.analyze_target("CC(=O)Nc1ccccc1", max_depth=5)
        assert tree.target == "CC(=O)Nc1ccccc1"
        assert len(tree.routes) > 0
        assert tree.num_alternatives == len(tree.routes)

    def test_mock_planner_respects_max_depth(self):
        """Routes don't exceed max_depth steps."""
        planner = RetrosynthesisPlanner(backend="mock")
        tree = planner.analyze_target("CC(=O)O", max_depth=3)
        for route in tree.routes:
            assert route.total_steps <= 3

    def test_score_route_valid_inputs(self):
        """Score route returns valid RouteScore."""
        mock_route = MagicMock(spec=SynthesisRoute)
        mock_route.total_steps = 3
        mock_route.confidence = 0.85
        mock_steps = []
        for _ in range(3):
            mock_step = MagicMock(spec=SynthesisStep)
            mock_step.reaction_type = "amide_coupling"
            mock_step.reactants = ["CC", "O"]
            mock_steps.append(mock_step)
        mock_route.steps = mock_steps
        
        score = score_route(mock_route)
        
        assert isinstance(score, RouteScore)
        assert hasattr(score, 'total_score')
        assert hasattr(score, 'step_count')
        assert hasattr(score, 'complexity_score')

    def test_score_route_confidence_bounds(self):
        """Confidence scores are 0.0-1.0."""
        mock_route = MagicMock(spec=SynthesisRoute)
        mock_route.total_steps = 2
        mock_route.confidence = 0.90
        mock_steps = []
        for _ in range(2):
            mock_step = MagicMock(spec=SynthesisStep)
            mock_step.reaction_type = "reduction"
            mock_step.reactants = ["CCO"]
            mock_steps.append(mock_step)
        mock_route.steps = mock_steps
        
        score = score_route(mock_route)
        
        assert 0.0 <= score.total_score <= 100.0
        assert score.step_count >= 0

    def test_check_building_blocks_integration(self):
        """Building block checker integrates with procurement."""
        smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
        
        availability = check_building_blocks(smiles_list)
        
        assert isinstance(availability, list)
        assert len(availability) == len(smiles_list)
        for result in availability:
            assert hasattr(result, 'smiles')
            assert hasattr(result, 'available')
            assert hasattr(result, 'vendors')

    def test_pluggable_backend_mock(self):
        """Backend='mock' works."""
        planner = RetrosynthesisPlanner(backend="mock")
        assert planner.backend == "mock"
        
        tree = planner.analyze_target("CC", max_depth=2)
        assert isinstance(tree, SynthesisTree)

    def test_pluggable_backend_unknown_raises(self):
        """Unknown backend raises ValueError."""
        planner = RetrosynthesisPlanner(backend="nonexistent")
        with pytest.raises(ValueError, match="Unknown backend"):
            planner.analyze_target("CCO", max_depth=3)

    def test_reaction_types_valid(self):
        """Reaction types are from valid taxonomy."""
        planner = RetrosynthesisPlanner(backend="mock")
        tree = planner.analyze_target("CC(=O)Nc1ccccc1", max_depth=3)
        
        valid_reaction_types = {
            "amide_coupling", "ester_hydrolysis", "friedel_crafts",
            "nucleophilic_substitution", "oxidation", "reduction",
            "alkylation", "acylation", "condensation"
        }
        
        for route in tree.routes:
            for step in route.steps:
                if hasattr(step, 'reaction_type'):
                    assert step.reaction_type in valid_reaction_types
