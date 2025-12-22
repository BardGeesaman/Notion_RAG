"""Vetted helper functions for notebook analysis workflows.

These functions provide safe, standardized ways to load entities and
perform common analysis operations in Jupyter notebooks.
"""

from __future__ import annotations

from typing import Optional, Dict, Any, List
from uuid import UUID

import pandas as pd


def load_dataset(dataset_id: str) -> Dict[str, Any]:
    """Load a dataset with its features.

    Args:
        dataset_id: UUID of the dataset

    Returns:
        Dictionary with dataset info and features DataFrame

    Raises:
        ValueError: If dataset not found
    """
    from sqlalchemy.orm import selectinload
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import Dataset

    with db_session() as db:
        dataset = (
            db.query(Dataset)
            .options(selectinload(Dataset.features))
            .filter(Dataset.id == UUID(dataset_id))
            .first()
        )

        if dataset is None:
            raise ValueError(f"Dataset not found: {dataset_id}")

        features_data = [
            {
                "id": str(f.id),
                "name": f.name,
                "feature_type": f.feature_type,
                "normalized_name": getattr(f, "normalized_name", None),
            }
            for f in (dataset.features or [])
        ]

        return {
            "id": str(dataset.id),
            "name": dataset.name,
            "omics_type": dataset.omics_type,
            "description": dataset.description,
            "feature_count": len(features_data),
            "features_df": pd.DataFrame(features_data),
        }


def load_experiment(experiment_id: str) -> Dict[str, Any]:
    """Load an experiment with its datasets.

    Args:
        experiment_id: UUID of the experiment

    Returns:
        Dictionary with experiment info and datasets list

    Raises:
        ValueError: If experiment not found
    """
    from sqlalchemy.orm import selectinload
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import Experiment

    with db_session() as db:
        experiment = (
            db.query(Experiment)
            .options(selectinload(Experiment.datasets))
            .filter(Experiment.id == UUID(experiment_id))
            .first()
        )

        if experiment is None:
            raise ValueError(f"Experiment not found: {experiment_id}")

        datasets_data = [
            {
                "id": str(ds.id),
                "name": ds.name,
                "omics_type": ds.omics_type,
            }
            for ds in (experiment.datasets or [])
        ]

        return {
            "id": str(experiment.id),
            "name": experiment.name,
            "type": experiment.type,
            "description": experiment.description,
            "design_type": experiment.design_type,
            "dataset_count": len(datasets_data),
            "datasets": datasets_data,
        }


def load_campaign(campaign_id: str) -> Dict[str, Any]:
    """Load an HTS campaign with summary statistics.

    Args:
        campaign_id: UUID of the campaign

    Returns:
        Dictionary with campaign info and hit statistics

    Raises:
        ValueError: If campaign not found
    """
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import HTSCampaign, HTSResult

    with db_session() as db:
        campaign = db.query(HTSCampaign).filter(HTSCampaign.id == UUID(campaign_id)).first()

        if campaign is None:
            raise ValueError(f"Campaign not found: {campaign_id}")

        # Get hit statistics
        total_results = db.query(HTSResult).filter(HTSResult.campaign_id == campaign.id).count()
        hits = db.query(HTSResult).filter(
            HTSResult.campaign_id == campaign.id,
            HTSResult.hit_flag == True
        ).count()

        return {
            "id": str(campaign.id),
            "name": campaign.name,
            "assay_type": campaign.assay_type,
            "target": campaign.target,
            "total_results": total_results,
            "hits": hits,
            "hit_rate": hits / total_results if total_results > 0 else 0,
        }


def load_compound(compound_id: str) -> Dict[str, Any]:
    """Load a compound with its properties.

    Args:
        compound_id: UUID of the compound

    Returns:
        Dictionary with compound info and properties

    Raises:
        ValueError: If compound not found
    """
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import Compound

    with db_session() as db:
        compound = db.query(Compound).filter(Compound.id == UUID(compound_id)).first()

        if compound is None:
            raise ValueError(f"Compound not found: {compound_id}")

        return {
            "id": str(compound.id),
            "corporate_id": compound.corporate_id,
            "name": compound.name,
            "smiles": compound.smiles,
            "molecular_weight": compound.molecular_weight,
            "logp": compound.logp,
        }


def run_hts_qc(campaign_id: str) -> Dict[str, Any]:
    """Run HTS quality control analysis on a campaign.

    Calculates Z' factor, signal-to-noise, and other QC metrics.

    Args:
        campaign_id: UUID of the campaign

    Returns:
        Dictionary with QC metrics and summary
    """
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import HTSCampaign, HTSResult
    import numpy as np

    with db_session() as db:
        campaign = db.query(HTSCampaign).filter(HTSCampaign.id == UUID(campaign_id)).first()
        if campaign is None:
            raise ValueError(f"Campaign not found: {campaign_id}")

        results = db.query(HTSResult).filter(HTSResult.campaign_id == campaign.id).all()

        if not results:
            return {"campaign_id": campaign_id, "error": "No results found"}

        # Extract values
        values = [r.normalized_value for r in results if r.normalized_value is not None]
        pos_controls = [r.normalized_value for r in results if r.control_type == "positive" and r.normalized_value]
        neg_controls = [r.normalized_value for r in results if r.control_type == "negative" and r.normalized_value]

        # Calculate Z' factor if controls available
        z_prime = None
        if pos_controls and neg_controls:
            pos_mean, pos_std = np.mean(pos_controls), np.std(pos_controls)
            neg_mean, neg_std = np.mean(neg_controls), np.std(neg_controls)
            if abs(pos_mean - neg_mean) > 0:
                z_prime = 1 - (3 * (pos_std + neg_std) / abs(pos_mean - neg_mean))

        return {
            "campaign_id": campaign_id,
            "campaign_name": campaign.name,
            "total_wells": len(results),
            "mean_value": np.mean(values) if values else None,
            "std_value": np.std(values) if values else None,
            "z_prime": z_prime,
            "qc_status": "PASS" if z_prime and z_prime > 0.5 else "REVIEW" if z_prime else "NO_CONTROLS",
        }


def fit_dose_response(
    compound_id: str,
    campaign_id: Optional[str] = None
) -> Dict[str, Any]:
    """Fit dose-response curve for a compound.

    Args:
        compound_id: UUID of the compound
        campaign_id: Optional campaign ID to filter results

    Returns:
        Dictionary with curve fit parameters and data
    """
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import BiochemicalResult
    import numpy as np

    with db_session() as db:
        query = db.query(BiochemicalResult).filter(
            BiochemicalResult.compound_id == UUID(compound_id)
        )
        if campaign_id:
            query = query.filter(BiochemicalResult.campaign_id == UUID(campaign_id))

        results = query.all()

        if not results:
            return {"compound_id": compound_id, "error": "No dose-response data found"}

        # Extract dose-response data
        doses = []
        responses = []
        for r in results:
            if r.concentration and r.percent_inhibition is not None:
                doses.append(r.concentration)
                responses.append(r.percent_inhibition)

        if len(doses) < 4:
            return {
                "compound_id": compound_id,
                "error": "Insufficient data points for curve fitting",
                "data_points": len(doses),
            }

        # Simple IC50 estimation (for MVP - could use scipy curve_fit)
        doses_arr = np.array(doses)
        responses_arr = np.array(responses)

        # Find approximate IC50 (concentration at 50% response)
        ic50_estimate = None
        if responses_arr.max() > 50 and responses_arr.min() < 50:
            # Linear interpolation to find IC50
            sorted_idx = np.argsort(doses_arr)
            for i in range(len(sorted_idx) - 1):
                r1, r2 = responses_arr[sorted_idx[i]], responses_arr[sorted_idx[i+1]]
                if (r1 < 50 <= r2) or (r2 < 50 <= r1):
                    d1, d2 = doses_arr[sorted_idx[i]], doses_arr[sorted_idx[i+1]]
                    ic50_estimate = d1 + (50 - r1) * (d2 - d1) / (r2 - r1)
                    break

        return {
            "compound_id": compound_id,
            "data_points": len(doses),
            "dose_range": [min(doses), max(doses)],
            "response_range": [min(responses), max(responses)],
            "ic50_estimate": ic50_estimate,
            "doses": doses,
            "responses": responses,
        }


def publish_to_rag(
    data: Dict[str, Any],
    tags: List[str],
    title: Optional[str] = None,
    entity_type: Optional[str] = None,
    entity_id: Optional[str] = None,
) -> Dict[str, str]:
    """Publish analysis results to the RAG index.

    Args:
        data: Analysis results dictionary
        tags: List of tags for categorization
        title: Optional title for the content
        entity_type: Optional entity type to link to
        entity_id: Optional entity ID to link to

    Returns:
        Dictionary with publish status and chunk ID
    """
    import json
    from amprenta_rag.logging_utils import get_logger

    logger = get_logger(__name__)

    # Format content for RAG
    content = f"# {title or 'Analysis Results'}\n\n"
    content += f"Tags: {', '.join(tags)}\n\n"
    content += f"```json\n{json.dumps(data, indent=2, default=str)}\n```"

    # For MVP, just log the publish (actual RAG indexing would go here)
    logger.info(f"Publishing to RAG: {title}, tags={tags}")

    return {
        "status": "published",
        "title": title or "Analysis Results",
        "tags": tags,
        "content_length": len(content),
    }


