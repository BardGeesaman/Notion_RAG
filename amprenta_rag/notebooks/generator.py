"""Jupyter notebook generators (nbformat v4)."""

from __future__ import annotations

from textwrap import dedent


def _md_cell(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source}


def _code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source,
    }


def generate_dataset_notebook(dataset_id: str) -> dict:
    """
    Generate an nbformat v4 notebook (as a dict) for exploring a dataset.

    Args:
        dataset_id: Dataset UUID (string)

    Returns:
        nbformat-compatible notebook dict (nbformat=4, nbformat_minor=5)
    """
    # Query dataset (name, omics_type, features) so the notebook has a good title.
    from uuid import UUID
    from sqlalchemy.orm import selectinload

    from amprenta_rag.database.models import Dataset
    from amprenta_rag.database.session import db_session

    dataset_uuid = UUID(dataset_id)
    with db_session() as db:
        dataset = (
            db.query(Dataset)
            .options(selectinload(Dataset.features))
            .filter(Dataset.id == dataset_uuid)
            .first()
        )

    if dataset is None:
        raise KeyError(f"Dataset not found: {dataset_id}")

    title_md = dedent(
        f"""\
        # Dataset Exploration: {dataset.name}

        **Dataset ID**: `{dataset_id}`
        **Omics Type**: `{dataset.omics_type}`
        **Feature Count**: `{len(dataset.features or [])}`
        """
    )

    imports_code = dedent(
        """\
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        sns.set_theme(style="whitegrid")
        """
    )

    db_code = dedent(
        """\
        from uuid import UUID
        from sqlalchemy.orm import selectinload

        from amprenta_rag.database.session import db_session
        from amprenta_rag.database.models import Dataset

        dataset_uuid = UUID(DATASET_ID)

        with db_session() as db:
            dataset = (
                db.query(Dataset)
                .options(selectinload(Dataset.features))
                .filter(Dataset.id == dataset_uuid)
                .first()
            )

        if dataset is None:
            raise ValueError(f"Dataset not found: {DATASET_ID}")

        print("Dataset:", dataset.name)
        print("Omics type:", dataset.omics_type)
        print("Feature count:", len(dataset.features or []))
        """
    )

    features_df_code = dedent(
        """\
        rows = []
        for f in (dataset.features or []):
            rows.append(
                {
                    "feature_id": str(f.id),
                    "name": f.name,
                    "feature_type": f.feature_type,
                    "normalized_name": getattr(f, "normalized_name", None),
                }
            )

        df = pd.DataFrame(rows)
        df
        """
    )

    describe_code = dedent(
        """\
        df.head()
        """
    ) + "\n" + dedent(
        """\
        df.describe(include="all")
        """
    )

    viz_code = dedent(
        """\
        if not df.empty and "feature_type" in df.columns:
            plt.figure(figsize=(8, 4))
            ax = sns.countplot(data=df, x="feature_type", order=df["feature_type"].value_counts().index)
            ax.set_title("Feature counts by type")
            ax.set_xlabel("Feature type")
            ax.set_ylabel("Count")
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            plt.show()
        else:
            print("No features found to visualize.")
        """
    )

    nb = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "name": "python3",
                "display_name": "Python 3",
                "language": "python",
            },
            "language_info": {
                "name": "python",
            },
        },
        "cells": [
            _md_cell(title_md),
            _code_cell(dedent(f'DATASET_ID = "{dataset_id}"\n')),
            _code_cell(imports_code),
            _code_cell(db_code),
            _code_cell(features_df_code),
            _code_cell(describe_code),
            _code_cell(viz_code),
        ],
    }

    return nb


def generate_context_cell(context) -> str:
    """Generate Python code for a notebook cell that sets analysis context."""
    from amprenta_rag.notebook.context import generate_context_cell as _generate_context_cell

    return _generate_context_cell(context)


def generate_hts_qc_cell(campaign_id: str) -> str:
    """Generate Python code for an HTS QC analysis cell."""
    return dedent(
        f"""\
        # HTS QC
        from amprenta_rag.notebook import run_hts_qc

        qc = run_hts_qc("{campaign_id}")
        qc
        """
    )


def generate_dose_response_cell(compound_id: str) -> str:
    """Generate Python code for a dose-response analysis cell."""
    return dedent(
        f"""\
        # Dose-response fit
        from amprenta_rag.notebook import fit_dose_response

        dr = fit_dose_response("{compound_id}")
        dr
        """
    )


def generate_publish_cell() -> str:
    """Generate Python code for publishing results to the RAG index."""
    return dedent(
        """\
        # Publish analysis to RAG (edit placeholders as needed)
        from amprenta_rag.notebook import publish_to_rag

        RESULTS = {}  # populate with your analysis output

        publish_result = publish_to_rag(
            data=RESULTS,
            tags=["analysis"],
            title="Analysis Results",
            entity_type=None,
            entity_id=None,
        )
        publish_result
        """
    )

