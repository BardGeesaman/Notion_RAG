"""
RAG text generation from Postgres data.

This module generates rich text summaries from Postgres entities for embedding
in the RAG system. Replaces Notion API-based text generation with efficient
Postgres queries using SQLAlchemy relationships.

Key Benefits:
- 10-100x faster than Notion API calls
- No rate limits
- Richer aggregations and statistics
- Offline-capable
- Consistent formatting
"""

from __future__ import annotations

from typing import Optional
from uuid import UUID

from sqlalchemy.orm import Session, selectinload

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Program, Experiment, Dataset, Feature, Signature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def build_program_text(db: Session, program_id: UUID) -> str:
    """
    Build comprehensive text summary for a program from Postgres.
    
    Args:
        db: Database session
        program_id: Program UUID
    
    Returns:
        Formatted text summary for RAG embedding
    """
    # Load program with all relationships
    program = db.query(Program).options(
        selectinload(Program.experiments),
        selectinload(Program.datasets),
        selectinload(Program.signatures),
    ).filter(Program.id == program_id).first()
    
    if not program:
        logger.warning("[RAG][TEXT-BUILDER] Program %s not found", program_id)
        return ""
    
    # Build text sections
    sections = []
    
    # Header
    sections.append(f"# Program: {program.name}")
    sections.append("")
    
    # Description
    if program.description:
        sections.append(f"**Description**: {program.description}")
        sections.append("")
    
    # Disease context
    if program.disease:
        sections.append(f"**Disease Areas**: {', '.join(program.disease)}")
        sections.append("")
    
    # Experiments
    if program.experiments:
        sections.append(f"## Experiments ({len(program.experiments)})")
        for exp in program.experiments:
            sections.append(f"- {exp.name}")
            if exp.type:
                sections.append(f"  - Type: {exp.type}")
            if exp.matrix:
                sections.append(f"  - Matrix: {', '.join(exp.matrix)}")
        sections.append("")
    
    # Datasets
    if program.datasets:
        # Group by omics type
        by_omics = {}
        for ds in program.datasets:
            by_omics.setdefault(ds.omics_type, []).append(ds)
        
        sections.append(f"## Datasets ({len(program.datasets)})")
        for omics_type, datasets in sorted(by_omics.items()):
            sections.append(f"### {omics_type.title()} ({len(datasets)} datasets)")
            for ds in datasets:
                sections.append(f"- {ds.name}")
                if ds.description:
                    sections.append(f"  - {ds.description[:200]}...")
        sections.append("")
    
    # Signatures
    if program.signatures:
        sections.append(f"## Signatures ({len(program.signatures)})")
        for sig in program.signatures:
            sections.append(f"- {sig.name}")
            if sig.modalities:
                sections.append(f"  - Modalities: {', '.join(sig.modalities)}")
        sections.append("")
    
    # Metadata
    sections.append("## Metadata")
    sections.append(f"- Program ID: {program.id}")
    sections.append(f"- Created: {program.created_at.isoformat()}")
    sections.append(f"- Updated: {program.updated_at.isoformat()}")
    
    return "\n".join(sections)


def build_dataset_text(db: Session, dataset_id: UUID) -> str:
    """
    Build comprehensive text summary for a dataset from Postgres.
    
    Args:
        db: Database session
        dataset_id: Dataset UUID
    
    Returns:
        Formatted text summary for RAG embedding
    """
    # Load dataset with relationships
    dataset = db.query(Dataset).options(
        selectinload(Dataset.programs),
        selectinload(Dataset.experiments),
        selectinload(Dataset.features),
        selectinload(Dataset.signatures),
    ).filter(Dataset.id == dataset_id).first()
    
    if not dataset:
        logger.warning("[RAG][TEXT-BUILDER] Dataset %s not found", dataset_id)
        return ""
    
    sections = []
    
    # Header
    sections.append(f"# Dataset: {dataset.name}")
    sections.append("")
    sections.append(f"**Omics Type**: {dataset.omics_type}")
    sections.append("")
    
    # Description
    if dataset.description:
        sections.append(f"**Description**: {dataset.description}")
        sections.append("")
    
    # Scientific metadata
    if dataset.organism:
        sections.append(f"**Organism**: {', '.join(dataset.organism)}")
    if dataset.sample_type:
        sections.append(f"**Sample Type**: {', '.join(dataset.sample_type)}")
    if dataset.disease:
        sections.append(f"**Disease**: {', '.join(dataset.disease)}")
    if dataset.organism or dataset.sample_type or dataset.disease:
        sections.append("")
    
    # Methods
    if dataset.methods:
        sections.append(f"## Methods")
        sections.append(dataset.methods)
        sections.append("")
    
    # Summary
    if dataset.summary:
        sections.append(f"## Summary")
        sections.append(dataset.summary)
        sections.append("")
    
    # Results
    if dataset.results:
        sections.append(f"## Results")
        sections.append(dataset.results)
        sections.append("")
    
    # Conclusions
    if dataset.conclusions:
        sections.append(f"## Conclusions")
        sections.append(dataset.conclusions)
        sections.append("")
    
    # Features
    if dataset.features:
        # Group by feature type
        by_type = {}
        for feat in dataset.features:
            by_type.setdefault(feat.feature_type, []).append(feat)
        
        sections.append(f"## Features ({len(dataset.features)} total)")
        for feat_type, features in sorted(by_type.items()):
            sections.append(f"- {feat_type}: {len(features)} features")
            # Show first 10 as examples
            if len(features) <= 10:
                for f in features:
                    sections.append(f"  - {f.name}")
            else:
                for f in features[:10]:
                    sections.append(f"  - {f.name}")
                sections.append(f"  - ... and {len(features) - 10} more")
        sections.append("")
    
    # Linked programs
    if dataset.programs:
        sections.append(f"## Programs")
        for prog in dataset.programs:
            sections.append(f"- {prog.name}")
        sections.append("")
    
    # Linked experiments
    if dataset.experiments:
        sections.append(f"## Experiments")
        for exp in dataset.experiments:
            sections.append(f"- {exp.name}")
        sections.append("")
    
    # Signature matches
    if dataset.signatures:
        sections.append(f"## Matched Signatures ({len(dataset.signatures)})")
        for sig in dataset.signatures:
            sections.append(f"- {sig.name}")
        sections.append("")
    
    # Metadata
    sections.append("## Metadata")
    sections.append(f"- Dataset ID: {dataset.id}")
    sections.append(f"- Data Origin: {dataset.data_origin or 'Internal'}")
    sections.append(f"- Source Type: {dataset.dataset_source_type or 'Unknown'}")
    sections.append(f"- Created: {dataset.created_at.isoformat()}")
    sections.append(f"- Updated: {dataset.updated_at.isoformat()}")
    if dataset.file_paths:
        sections.append(f"- Files: {', '.join(dataset.file_paths)}")
    
    return "\n".join(sections)


def build_signature_text(db: Session, signature_id: UUID) -> str:
    """
    Build comprehensive text summary for a signature from Postgres.
    
    Args:
        db: Database session
        signature_id: Signature UUID
    
    Returns:
        Formatted text summary for RAG embedding
    """
    # Load signature with relationships
    signature = db.query(Signature).options(
        selectinload(Signature.components),
        selectinload(Signature.programs),
        selectinload(Signature.datasets),
    ).filter(Signature.id == signature_id).first()
    
    if not signature:
        logger.warning("[RAG][TEXT-BUILDER] Signature %s not found", signature_id)
        return ""
    
    sections = []
    
    # Header
    sections.append(f"# Signature: {signature.name}")
    sections.append("")
    
    # Description
    if signature.description:
        sections.append(f"**Description**: {signature.description}")
        sections.append("")
    
    # Metadata
    if signature.modalities:
        sections.append(f"**Modalities**: {', '.join(signature.modalities)}")
    if signature.short_id:
        sections.append(f"**Short ID**: {signature.short_id}")
    if signature.biomarker_role:
        sections.append(f"**Biomarker Role**: {', '.join(signature.biomarker_role)}")
    if signature.phenotype_axes:
        sections.append(f"**Phenotype Axes**: {', '.join(signature.phenotype_axes)}")
    if signature.data_ownership:
        sections.append(f"**Data Ownership**: {signature.data_ownership}")
    sections.append("")
    
    # Components
    if signature.components:
        # Group by feature type
        by_type = {}
        for comp in signature.components:
            by_type.setdefault(comp.feature_type, []).append(comp)
        
        sections.append(f"## Components ({len(signature.components)} total)")
        for feat_type, components in sorted(by_type.items()):
            sections.append(f"### {feat_type.title()} ({len(components)} components)")
            for comp in components:
                direction = comp.direction or "neutral"
                weight = comp.weight or 1.0
                sections.append(f"- {comp.feature_name} ({direction}, weight: {weight:.2f})")
        sections.append("")
    
    # Linked programs
    if signature.programs:
        sections.append(f"## Programs")
        for prog in signature.programs:
            sections.append(f"- {prog.name}")
        sections.append("")
    
    # Matched datasets
    if signature.datasets:
        sections.append(f"## Matched Datasets ({len(signature.datasets)})")
        for ds in signature.datasets:
            sections.append(f"- {ds.name} ({ds.omics_type})")
        sections.append("")
    
    # Metadata
    sections.append("## Metadata")
    sections.append(f"- Signature ID: {signature.id}")
    sections.append(f"- Created: {signature.created_at.isoformat()}")
    sections.append(f"- Updated: {signature.updated_at.isoformat()}")
    
    return "\n".join(sections)


def build_feature_text(db: Session, feature_id: UUID) -> str:
    """
    Build text summary for a feature from Postgres.
    
    Args:
        db: Database session
        feature_id: Feature UUID
    
    Returns:
        Formatted text summary for RAG embedding
    """
    # Load feature with relationships
    feature = db.query(Feature).options(
        selectinload(Feature.datasets),
        selectinload(Feature.signatures),
    ).filter(Feature.id == feature_id).first()
    
    if not feature:
        logger.warning("[RAG][TEXT-BUILDER] Feature %s not found", feature_id)
        return ""
    
    sections = []
    
    # Header
    sections.append(f"# Feature: {feature.name}")
    sections.append("")
    sections.append(f"**Type**: {feature.feature_type}")
    sections.append("")
    
    # Normalized name and aliases
    if feature.normalized_name and feature.normalized_name != feature.name:
        sections.append(f"**Normalized Name**: {feature.normalized_name}")
    if feature.aliases:
        sections.append(f"**Aliases**: {', '.join(feature.aliases)}")
    if feature.normalized_name or feature.aliases:
        sections.append("")
    
    # External IDs
    if feature.external_ids:
        sections.append("## External IDs")
        for key, value in feature.external_ids.items():
            sections.append(f"- {key}: {value}")
        sections.append("")
    
    # Datasets
    if feature.datasets:
        # Group by omics type
        by_omics = {}
        for ds in feature.datasets:
            by_omics.setdefault(ds.omics_type, []).append(ds)
        
        sections.append(f"## Datasets ({len(feature.datasets)} total)")
        for omics_type, datasets in sorted(by_omics.items()):
            sections.append(f"### {omics_type.title()} ({len(datasets)} datasets)")
            for ds in datasets[:5]:  # Show first 5
                sections.append(f"- {ds.name}")
            if len(datasets) > 5:
                sections.append(f"- ... and {len(datasets) - 5} more")
        sections.append("")
    
    # Signatures
    if feature.signatures:
        sections.append(f"## Signatures ({len(feature.signatures)})")
        for sig in feature.signatures:
            sections.append(f"- {sig.name}")
        sections.append("")
    
    # Metadata
    sections.append("## Metadata")
    sections.append(f"- Feature ID: {feature.id}")
    sections.append(f"- Created: {feature.created_at.isoformat()}")
    sections.append(f"- Updated: {feature.updated_at.isoformat()}")
    
    return "\n".join(sections)


def build_experiment_text(db: Session, experiment_id: UUID) -> str:
    """
    Build text summary for an experiment from Postgres.
    
    Args:
        db: Database session
        experiment_id: Experiment UUID
    
    Returns:
        Formatted text summary for RAG embedding
    """
    # Load experiment with relationships
    experiment = db.query(Experiment).options(
        selectinload(Experiment.programs),
        selectinload(Experiment.datasets),
    ).filter(Experiment.id == experiment_id).first()
    
    if not experiment:
        logger.warning("[RAG][TEXT-BUILDER] Experiment %s not found", experiment_id)
        return ""
    
    sections = []
    
    # Header
    sections.append(f"# Experiment: {experiment.name}")
    sections.append("")
    
    # Type and description
    if experiment.type:
        sections.append(f"**Type**: {experiment.type}")
    if experiment.description:
        sections.append(f"**Description**: {experiment.description}")
    sections.append("")
    
    # Scientific metadata
    if experiment.disease:
        sections.append(f"**Disease**: {', '.join(experiment.disease)}")
    if experiment.matrix:
        sections.append(f"**Matrix**: {', '.join(experiment.matrix)}")
    if experiment.model_systems:
        sections.append(f"**Model Systems**: {', '.join(experiment.model_systems)}")
    if experiment.targets:
        sections.append(f"**Targets**: {', '.join(experiment.targets)}")
    if experiment.modality:
        sections.append(f"**Modality**: {', '.join(experiment.modality)}")
    if experiment.stage:
        sections.append(f"**Stage**: {experiment.stage}")
    if experiment.biomarker_role:
        sections.append(f"**Biomarker Role**: {', '.join(experiment.biomarker_role)}")
    if experiment.treatment_arms:
        sections.append(f"**Treatment Arms**: {', '.join(experiment.treatment_arms)}")
    sections.append("")
    
    # Programs
    if experiment.programs:
        sections.append(f"## Programs")
        for prog in experiment.programs:
            sections.append(f"- {prog.name}")
        sections.append("")
    
    # Datasets
    if experiment.datasets:
        # Group by omics type
        by_omics = {}
        for ds in experiment.datasets:
            by_omics.setdefault(ds.omics_type, []).append(ds)
        
        sections.append(f"## Datasets ({len(experiment.datasets)} total)")
        for omics_type, datasets in sorted(by_omics.items()):
            sections.append(f"### {omics_type.title()} ({len(datasets)} datasets)")
            for ds in datasets:
                sections.append(f"- {ds.name}")
        sections.append("")
    
    # Metadata
    sections.append("## Metadata")
    sections.append(f"- Experiment ID: {experiment.id}")
    sections.append(f"- Created: {experiment.created_at.isoformat()}")
    sections.append(f"- Updated: {experiment.updated_at.isoformat()}")
    
    return "\n".join(sections)


# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

def build_text_for_entity(
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None,
) -> str:
    """
    Build text for any entity type.
    
    Args:
        entity_type: One of: program, experiment, dataset, feature, signature
        entity_id: Entity UUID
        db: Database session (optional)
    
    Returns:
        Formatted text summary
    """
    if db is None:
        db = next(get_db())
        close_db = True
    else:
        close_db = False
    
    try:
        entity_type_lower = entity_type.lower()
        
        if entity_type_lower == "program":
            return build_program_text(db, entity_id)
        elif entity_type_lower == "experiment":
            return build_experiment_text(db, entity_id)
        elif entity_type_lower == "dataset":
            return build_dataset_text(db, entity_id)
        elif entity_type_lower == "feature":
            return build_feature_text(db, entity_id)
        elif entity_type_lower == "signature":
            return build_signature_text(db, entity_id)
        else:
            logger.error("[RAG][TEXT-BUILDER] Unknown entity type: %s", entity_type)
            return ""
    
    finally:
        if close_db:
            db.close()

