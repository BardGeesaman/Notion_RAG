"""Retrosynthesis planning service for synthesis route analysis."""

from __future__ import annotations

import uuid
from dataclasses import dataclass
from typing import List


@dataclass
class SynthesisStep:
    reactants: List[str]       # SMILES of reactants
    product: str               # SMILES of product
    reaction_type: str         # e.g., "amide_coupling", "reduction"
    conditions: str            # Reaction conditions
    confidence: float          # 0.0-1.0


@dataclass
class SynthesisRoute:
    id: str                    # UUID for caching
    steps: List[SynthesisStep]
    total_steps: int
    confidence: float


@dataclass
class SynthesisTree:
    target: str                # Target SMILES
    routes: List[SynthesisRoute]
    analysis_time_ms: int
    num_alternatives: int


@dataclass
class RouteScore:
    total_score: float         # 0-100
    step_count: int
    complexity_score: float    # Lower = simpler
    availability_score: float  # Higher = more available
    estimated_cost: float      # USD estimate


@dataclass
class BuildingBlockResult:
    smiles: str
    available: bool
    vendors: List[dict]


class RetrosynthesisPlanner:
    """Retrosynthesis route planner with pluggable backend."""
    
    def __init__(self, backend: str = "mock"):
        self.backend = backend  # "mock" | "ibm_rxn" | "aizynthfinder"
    
    def analyze_target(self, smiles: str, max_depth: int = 5) -> SynthesisTree:
        """Analyze target molecule and return synthesis routes."""
        if self.backend == "mock":
            return self._mock_analysis(smiles, max_depth)
        # Future: elif self.backend == "ibm_rxn": ...
        raise ValueError(f"Unknown backend: {self.backend}")
    
    def _mock_analysis(self, smiles: str, max_depth: int) -> SynthesisTree:
        """Return mock synthesis tree for demo/testing."""
        import time
        start_time = time.time()
        
        # Generate 2-3 alternative routes
        routes = []
        
        # Route 1: Traditional coupling approach
        route1_steps = [
            SynthesisStep(
                reactants=["CC(C)OC(=O)c1ccc(N)cc1", "CC(=O)Cl"],
                product="CC(C)OC(=O)c1ccc(NC(=O)C)cc1",
                reaction_type="amide_coupling",
                conditions="EDC, HOBt, DMF, rt, 2h",
                confidence=0.85
            ),
            SynthesisStep(
                reactants=["CC(C)OC(=O)c1ccc(NC(=O)C)cc1"],
                product=smiles,
                reaction_type="ester_hydrolysis",
                conditions="LiOH, THF/H2O, rt, 4h",
                confidence=0.90
            )
        ]
        
        route1 = SynthesisRoute(
            id=str(uuid.uuid4()),
            steps=route1_steps,
            total_steps=len(route1_steps),
            confidence=0.87
        )
        routes.append(route1)
        
        # Route 2: Alternative approach
        route2_steps = [
            SynthesisStep(
                reactants=["Nc1ccc(C(=O)O)cc1", "CC(=O)Cl"],
                product="CC(=O)Nc1ccc(C(=O)O)cc1",
                reaction_type="amide_coupling",
                conditions="PyBOP, DIPEA, DCM, rt, 1h",
                confidence=0.92
            )
        ]
        
        route2 = SynthesisRoute(
            id=str(uuid.uuid4()),
            steps=route2_steps,
            total_steps=len(route2_steps),
            confidence=0.92
        )
        routes.append(route2)
        
        # Route 3: Multi-step approach
        route3_steps = [
            SynthesisStep(
                reactants=["CC(C)OC(=O)c1ccc(Br)cc1", "CC(=O)NHNH2"],
                product="CC(C)OC(=O)c1ccc(NHNC(=O)C)cc1",
                reaction_type="buchwald_hartwig",
                conditions="Pd(OAc)2, BINAP, Cs2CO3, toluene, 110°C, 12h",
                confidence=0.75
            ),
            SynthesisStep(
                reactants=["CC(C)OC(=O)c1ccc(NHNC(=O)C)cc1"],
                product="CC(C)OC(=O)c1ccc(NC(=O)C)cc1",
                reaction_type="reduction",
                conditions="Zn, AcOH, rt, 2h",
                confidence=0.80
            ),
            SynthesisStep(
                reactants=["CC(C)OC(=O)c1ccc(NC(=O)C)cc1"],
                product=smiles,
                reaction_type="ester_hydrolysis",
                conditions="NaOH, MeOH, reflux, 3h",
                confidence=0.88
            )
        ]
        
        route3 = SynthesisRoute(
            id=str(uuid.uuid4()),
            steps=route3_steps,
            total_steps=len(route3_steps),
            confidence=0.81
        )
        routes.append(route3)
        
        analysis_time = int((time.time() - start_time) * 1000)
        
        return SynthesisTree(
            target=smiles,
            routes=routes,
            analysis_time_ms=analysis_time,
            num_alternatives=len(routes)
        )


def score_route(route: SynthesisRoute) -> RouteScore:
    """Score a synthesis route based on complexity and availability."""
    # Mock scoring algorithm
    step_count = route.total_steps
    
    # Complexity: fewer steps = lower complexity
    complexity_score = max(0, 100 - (step_count * 15))
    
    # Availability: based on reaction types (mock)
    common_reactions = {"amide_coupling", "ester_hydrolysis", "reduction"}
    availability_score = sum(
        80 if step.reaction_type in common_reactions else 60
        for step in route.steps
    ) / len(route.steps)
    
    # Cost: rough estimate based on steps
    estimated_cost = step_count * 150 + sum(
        len(step.reactants) * 50 for step in route.steps
    )
    
    # Total score: weighted combination
    total_score = (
        complexity_score * 0.3 +
        availability_score * 0.4 +
        route.confidence * 100 * 0.3
    )
    
    return RouteScore(
        total_score=total_score,
        step_count=step_count,
        complexity_score=complexity_score,
        availability_score=availability_score,
        estimated_cost=estimated_cost
    )


def check_building_blocks(smiles_list: List[str]) -> List[BuildingBlockResult]:
    """Check building block availability via procurement service."""
    results = []
    
    for smiles in smiles_list:
        # Mock availability check
        # In reality, would integrate with procurement service
        is_available = len(smiles) < 30  # Simple heuristic for mock
        
        vendors = []
        if is_available:
            vendors = [
                {"name": "Sigma-Aldrich", "catalog_id": f"SA{hash(smiles) % 100000}", "price": "$45-150"},
                {"name": "TCI", "catalog_id": f"T{hash(smiles) % 10000}", "price": "$35-120"},
            ]
        
        results.append(BuildingBlockResult(
            smiles=smiles,
            available=is_available,
            vendors=vendors
        ))
    
    return results
