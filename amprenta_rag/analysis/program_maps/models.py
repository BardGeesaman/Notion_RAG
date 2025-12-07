"""
Data models for program-signature mapping.

Defines data classes for representing program-signature relationships,
omics coverage, and convergence metrics.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class ProgramSignatureScore:
    """
    Represents a signature score for a program.

    Attributes:
        program_id: Notion page ID of program
        signature_id: Notion page ID of signature
        program_name: Name of program
        signature_name: Name of signature
        overall_score: Overall signature score (0-1)
        score_by_omics: Scores by omics type
        matching_datasets: List of dataset IDs that match this signature
        coverage_fraction: Fraction of program datasets that match
    """

    program_id: str
    signature_id: str
    program_name: str
    signature_name: str
    overall_score: float
    score_by_omics: Dict[str, float] = field(default_factory=dict)
    matching_datasets: List[str] = field(default_factory=list)
    coverage_fraction: float = 0.0


@dataclass
class ProgramOmicsCoverage:
    """
    Represents omics coverage for a program.

    Attributes:
        program_id: Notion page ID of program
        program_name: Name of program
        total_datasets: Total number of datasets in program
        datasets_by_omics: Number of datasets by omics type
        features_by_omics: Number of unique features by omics type
        coverage_summary: Summary of omics coverage
    """

    program_id: str
    program_name: str
    total_datasets: int
    datasets_by_omics: Dict[str, int] = field(default_factory=dict)
    features_by_omics: Dict[str, int] = field(default_factory=dict)
    coverage_summary: str = ""


@dataclass
class ProgramSignatureMap:
    """
    Represents a complete program-signature mapping.

    Attributes:
        program_id: Notion page ID of program
        program_name: Name of program
        signature_scores: List of ProgramSignatureScore objects
        omics_coverage: ProgramOmicsCoverage object
        top_signatures: Top N signatures by score
        convergence_indicators: Cross-omics convergence metrics
    """

    program_id: str
    program_name: str
    signature_scores: List[ProgramSignatureScore] = field(default_factory=list)
    omics_coverage: Optional[ProgramOmicsCoverage] = None
    top_signatures: List[ProgramSignatureScore] = field(default_factory=list)
    convergence_indicators: Dict[str, float] = field(default_factory=dict)
