from __future__ import annotations

"""
Utility script to seed Postgres with a demo dataset and features
containing differential expression statistics for visualization tests.
"""

from uuid import uuid4

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset, Feature, Program


def create_test_visualization_data() -> None:
    db_gen = get_db()
    db = next(db_gen)

    dataset_name = "Test DE Analysis - ALS vs Control"
    program_name = "Test Visualization Program"

    try:
        # Ensure program exists
        program = db.query(Program).filter(Program.name == program_name).first()
        if not program:
            program = Program(
                id=uuid4(),
                name=program_name,
                description="Demo data for visualization testing",
                disease=["ALS"],
            )
            db.add(program)

        # Prevent duplicates if rerun
        existing = db.query(Dataset).filter(Dataset.name == dataset_name).first()
        if existing:
            print(f"Dataset already exists: {existing.id} ({existing.name})")
            return

        dataset = Dataset(
            id=uuid4(),
            name=dataset_name,
            omics_type="transcriptomics",
            description="Simulated differential expression data",
            disease=["ALS"],
            sample_type=["CSF"],
        )
        db.add(dataset)
        dataset.programs.append(program)

        # Define 50 genes with varying stats
        significant_up = [
            ("TP53", 2.5, 0.0001),
            ("TNF", 1.9, 0.0010),
            ("BRCA1", 1.7, 0.0050),
            ("APP", 2.1, 0.0008),
            ("MAPT", 1.6, 0.0060),
            ("SOD1", 2.2, 0.0006),
            ("TARDBP", 1.8, 0.0090),
            ("FUS", 2.0, 0.0040),
            ("C9ORF72", 1.55, 0.0070),
            ("GRN", 1.65, 0.0050),
        ]

        significant_down = [
            ("IL6", -2.1, 0.0009),
            ("CXCL8", -1.9, 0.0020),
            ("HSPA1A", -1.7, 0.0040),
            ("HSP90AA1", -1.6, 0.0050),
            ("UBQLN2", -1.8, 0.0030),
            ("SQSTM1", -2.0, 0.0007),
            ("OPTN", -1.55, 0.0060),
            ("VCP", -1.75, 0.0045),
            ("ATXN2", -1.6, 0.0080),
            ("CHMP2B", -1.9, 0.0015),
        ]

        moderate = [
            ("GATA3", 1.1, 0.020),
            ("STAT3", 0.9, 0.030),
            ("NFKB1", 1.3, 0.015),
            ("JAK2", 0.8, 0.040),
            ("PIK3CA", 1.0, 0.025),
            ("AKT1", 1.4, 0.018),
            ("EGFR", 0.7, 0.032),
            ("VEGFA", 1.2, 0.022),
            ("MMP9", 1.25, 0.017),
            ("ITGAM", 0.95, 0.028),
            ("CXCR4", -1.1, 0.021),
            ("CCR5", -1.2, 0.019),
            ("SLC1A2", -0.9, 0.033),
            ("GLUL", -1.3, 0.016),
            ("ENO1", -0.8, 0.035),
            ("LDHA", -1.25, 0.027),
            ("PKM", -1.4, 0.014),
            ("HIF1A", -1.0, 0.026),
            ("SIRT1", 1.05, 0.029),
            ("MTOR", 0.85, 0.034),
        ]

        not_significant = [
            ("ACTB", 0.2, 0.20),
            ("GAPDH", -0.1, 0.50),
            ("RPLP0", 0.3, 0.15),
            ("B2M", -0.2, 0.40),
            ("PPIA", 0.1, 0.60),
            ("EEF1A1", -0.3, 0.25),
            ("HPRT1", 0.05, 0.80),
            ("TUBB", -0.15, 0.55),
            ("RPL13A", 0.25, 0.18),
            ("YWHAZ", -0.05, 0.70),
        ]

        gene_stats = significant_up + significant_down + moderate + not_significant

        for gene_name, log2fc, pval in gene_stats:
            feature = Feature(
                id=uuid4(),
                name=gene_name,
                feature_type="gene",
                external_ids={
                    "log2FC": log2fc,
                    "log2FoldChange": log2fc,
                    "fold_change": log2fc,
                    "pvalue": pval,
                    "p_val": pval,
                    "padj": min(pval * 50, 1.0),
                    "stats": {"log2FC": log2fc, "pvalue": pval},
                },
            )
            db.add(feature)
            dataset.features.append(feature)

        db.commit()
        print(f"✅ Created test dataset: {dataset.id}")
        print(f"   Name: {dataset.name}")
        print(f"   Features: {len(gene_stats)}")

    finally:
        try:
            db.close()
        finally:
            # Exhaust generator to trigger cleanup in get_db
            try:
                next(db_gen, None)
            except Exception:
                pass


if __name__ == "__main__":
    create_test_visualization_data()

