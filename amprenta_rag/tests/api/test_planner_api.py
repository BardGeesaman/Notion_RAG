"""Tests for experiment planner API endpoints."""

from __future__ import annotations

import asyncio
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestPowerEndpoint:
    """Tests for POST /api/v1/planner/power endpoint."""

    @patch("amprenta_rag.api.routers.planner.calculate_sample_size")
    def test_calculate_power_success(self, mock_calc):
        """Test successful power calculation."""
        mock_calc.return_value = 64
        
        response = client.post(
            "/api/v1/planner/power",
            json={
                "effect_size": 0.5,
                "alpha": 0.05,
                "power": 0.80,
                "test_type": "t-test",
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["n_per_group"] == 64


class TestEffectSizeEndpoint:
    """Tests for POST /api/v1/planner/effect-size endpoint."""

    @patch("amprenta_rag.api.routers.planner.estimate_effect_size_from_data")
    def test_estimate_effect_size_success(self, mock_estimate):
        """Test successful effect size estimation."""
        mock_estimate.return_value = 0.75
        
        response = client.post(
            "/api/v1/planner/effect-size",
            json={
                "group1": [10.0, 12.0, 11.0],
                "group2": [8.0, 9.0, 7.0],
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["effect_size"] == 0.75


class TestPlatesEndpoint:
    """Tests for POST /api/v1/planner/plates endpoint."""

    @patch("amprenta_rag.api.routers.planner.calculate_plate_layout")
    def test_calculate_plates_success(self, mock_calc):
        """Test successful plate layout calculation."""
        mock_calc.return_value = {
            "plates_needed": 2,
            "wells_used": 100,
            "empty_wells": 92,
        }
        
        response = client.post(
            "/api/v1/planner/plates",
            json={"n": 100, "plate_format": 96},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["plates_needed"] == 2
        assert data["wells_used"] == 100


class TestCostEndpoint:
    """Tests for POST /api/v1/planner/cost endpoint."""

    @patch("amprenta_rag.api.routers.planner.estimate_experiment_cost")
    def test_estimate_cost_success(self, mock_estimate):
        """Test successful cost estimation."""
        mock_estimate.return_value = {
            "sample_cost": 1000.0,
            "overhead": 100.0,
            "total": 1100.0,
        }
        
        response = client.post(
            "/api/v1/planner/cost",
            json={
                "n": 100,
                "cost_per_sample": 10.0,
                "overhead_pct": 0.1,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 1100.0

    @patch("amprenta_rag.api.routers.planner.estimate_experiment_cost")
    def test_estimate_cost_custom_overhead(self, mock_estimate):
        """Test cost estimation with custom overhead."""
        mock_estimate.return_value = {
            "sample_cost": 500.0,
            "overhead": 75.0,
            "total": 575.0,
        }
        
        response = client.post(
            "/api/v1/planner/cost",
            json={
                "n": 50,
                "cost_per_sample": 10.0,
                "overhead_pct": 0.15,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["overhead"] == 75.0


class TestAsyncLLMPlanning:
    """Test async execution of LLM-based planning endpoints."""

    @pytest.mark.asyncio
    async def test_plan_async(self):
        """Test plan endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.planner._sync_create_plan') as mock_plan:
            from amprenta_rag.api.schemas import PlanResult, PlanStep
            
            # Mock the planning function
            mock_result = PlanResult(
                plan_id="plan_123",
                title="ALS Drug Discovery Protocol",
                objective="Identify potential ALS therapeutic compounds",
                steps=[
                    PlanStep(
                        step_number=1,
                        description="Literature review and target identification",
                        duration="2 weeks",
                        resources=["PubMed access", "Research analyst"],
                        dependencies=[]
                    ),
                    PlanStep(
                        step_number=2,
                        description="High-throughput screening setup",
                        duration="1 week",
                        resources=["384-well plates", "Compound library"],
                        dependencies=[1]
                    )
                ],
                estimated_duration="6 months",
                estimated_cost="$50,000",
                risks=["Compound availability", "Assay validation"],
                processing_time_seconds=2.1,
                cached=False
            )
            mock_plan.return_value = mock_result
            
            from amprenta_rag.api.routers.planner import plan_endpoint
            from amprenta_rag.api.schemas import PlanRequest, ScoringContextRequest
            
            request = PlanRequest(
                goal="Develop a screening protocol for ALS therapeutic compounds",
                context=ScoringContextRequest(
                    diseases=["ALS"],
                    targets=["SOD1", "TDP-43"],
                    species=["human"],
                    assay_types=["cell viability", "neuroprotection"],
                    min_sample_size=100,
                ),
                constraints=["Budget: $100K", "Timeline: 6 months"]
            )
            
            result = await plan_endpoint(request)
            
            # Verify async execution and result
            assert result.plan_id == "plan_123"
            assert result.title == "ALS Drug Discovery Protocol"
            assert len(result.steps) == 2
            assert result.steps[0].step_number == 1
            assert result.steps[1].dependencies == [1]
            assert result.estimated_duration == "6 months"
            assert "Compound availability" in result.risks
            mock_plan.assert_called_once()

    @pytest.mark.asyncio
    async def test_critique_async(self):
        """Test critique endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.planner._sync_critique_plan') as mock_critique:
            from amprenta_rag.api.schemas import CritiqueResult
            
            # Mock the critique function
            mock_result = CritiqueResult(
                overall_score=7.5,
                strengths=[
                    "Well-defined objectives",
                    "Appropriate statistical power"
                ],
                weaknesses=[
                    "Limited control conditions",
                    "Insufficient sample size for subgroup analysis"
                ],
                recommendations=[
                    "Add negative control groups",
                    "Increase sample size by 20%",
                    "Include dose-response curves"
                ],
                feasibility_score=8.0,
                processing_time_seconds=1.8,
                cached=False
            )
            mock_critique.return_value = mock_result
            
            from amprenta_rag.api.routers.planner import critique_endpoint
            from amprenta_rag.api.schemas import CritiqueRequest, ScoringContextRequest
            
            request = CritiqueRequest(
                plan="1. Screen 1000 compounds\n2. Test top 10 hits\n3. Validate in animal models",
                criteria=["scientific rigor", "feasibility", "cost-effectiveness"],
                context=ScoringContextRequest(
                    diseases=["ALS"],
                    targets=["SOD1"],
                    species=["human", "mouse"],
                    assay_types=["cell viability"],
                    min_sample_size=50,
                )
            )
            
            result = await critique_endpoint(request)
            
            # Verify async execution and result
            assert result.overall_score == 7.5
            assert len(result.strengths) == 2
            assert len(result.weaknesses) == 2
            assert len(result.recommendations) == 3
            assert result.feasibility_score == 8.0
            assert "Well-defined objectives" in result.strengths
            assert "Add negative control groups" in result.recommendations
            mock_critique.assert_called_once()

    @pytest.mark.asyncio
    async def test_refine_async(self):
        """Test refine endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.planner._sync_refine_plan') as mock_refine:
            from amprenta_rag.api.schemas import RefinementResult, PlanResult, PlanStep
            
            # Mock the refinement function
            refined_plan = PlanResult(
                plan_id="plan_123_v2",
                title="Improved ALS Drug Discovery Protocol",
                objective="Identify potential ALS therapeutic compounds with enhanced controls",
                steps=[
                    PlanStep(
                        step_number=1,
                        description="Literature review and target identification",
                        duration="2 weeks",
                        resources=["PubMed access", "Research analyst"],
                        dependencies=[]
                    ),
                    PlanStep(
                        step_number=2,
                        description="High-throughput screening with controls",
                        duration="2 weeks",
                        resources=["384-well plates", "Compound library", "Control compounds"],
                        dependencies=[1]
                    ),
                    PlanStep(
                        step_number=3,
                        description="Dose-response validation",
                        duration="1 week",
                        resources=["96-well plates", "Selected hits"],
                        dependencies=[2]
                    )
                ],
                estimated_duration="7 months",
                estimated_cost="$65,000",
                risks=["Compound availability"],
                processing_time_seconds=2.3,
                cached=False
            )
            
            mock_result = RefinementResult(
                refined_plan=refined_plan,
                changes_made=[
                    "Added negative control groups",
                    "Included dose-response validation step",
                    "Increased timeline to accommodate additional controls"
                ],
                rationale="Enhanced experimental rigor through proper controls and validation",
                processing_time_seconds=2.3,
                cached=False
            )
            mock_refine.return_value = mock_result
            
            from amprenta_rag.api.routers.planner import refine_endpoint
            from amprenta_rag.api.schemas import RefineRequest
            
            request = RefineRequest(
                original_plan="1. Screen 1000 compounds\n2. Test top 10 hits",
                critique="Needs better controls and validation steps",
                additional_requirements="Include dose-response analysis"
            )
            
            result = await refine_endpoint(request)
            
            # Verify async execution and result
            assert result.refined_plan.plan_id == "plan_123_v2"
            assert len(result.refined_plan.steps) == 3
            assert len(result.changes_made) == 3
            assert "Added negative control groups" in result.changes_made
            assert "Enhanced experimental rigor" in result.rationale
            assert result.refined_plan.estimated_cost == "$65,000"
            mock_refine.assert_called_once()

    @pytest.mark.asyncio
    async def test_execute_async(self):
        """Test execute endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.planner._sync_execute_plan') as mock_execute:
            from amprenta_rag.api.schemas import ExecutionGuidance, PlanStep
            
            # Mock the execution guidance function
            mock_result = ExecutionGuidance(
                current_step=PlanStep(
                    step_number=2,
                    description="High-throughput screening with controls",
                    duration="2 weeks",
                    resources=["384-well plates", "Compound library"],
                    dependencies=[1]
                ),
                next_steps=[
                    PlanStep(
                        step_number=3,
                        description="Hit validation",
                        duration="1 week",
                        resources=["96-well plates"],
                        dependencies=[2]
                    ),
                    PlanStep(
                        step_number=4,
                        description="Secondary assays",
                        duration="2 weeks",
                        resources=["Specialized equipment"],
                        dependencies=[3]
                    )
                ],
                recommendations=[
                    "Prepare validation plates in advance",
                    "Order specialized reagents for secondary assays",
                    "Schedule equipment time for next week"
                ],
                potential_issues=[
                    "Compound solubility problems",
                    "Equipment availability conflicts"
                ],
                troubleshooting=[
                    "Use DMSO as alternative solvent",
                    "Book backup equipment slots"
                ],
                processing_time_seconds=1.5,
                cached=False
            )
            mock_execute.return_value = mock_result
            
            from amprenta_rag.api.routers.planner import execute_endpoint
            from amprenta_rag.api.schemas import ExecuteRequest
            
            request = ExecuteRequest(
                plan="1. Literature review\n2. HTS screening\n3. Hit validation\n4. Secondary assays",
                current_step=2,
                issues=["Some compounds not dissolving properly"]
            )
            
            result = await execute_endpoint(request)
            
            # Verify async execution and result
            assert result.current_step.step_number == 2
            assert result.current_step.description == "High-throughput screening with controls"
            assert len(result.next_steps) == 2
            assert result.next_steps[0].step_number == 3
            assert len(result.recommendations) == 3
            assert "Prepare validation plates in advance" in result.recommendations
            assert len(result.potential_issues) == 2
            assert "Compound solubility problems" in result.potential_issues
            assert "Use DMSO as alternative solvent" in result.troubleshooting
            mock_execute.assert_called_once()

