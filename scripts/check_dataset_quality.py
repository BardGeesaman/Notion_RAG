from __future__ import annotations

from amprenta_rag.analysis.quality_metrics import compute_quality_score
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset


def main():
    db_gen = get_db()
    db = next(db_gen)
    try:
        datasets = db.query(Dataset).all()
        if not datasets:
            print("No datasets found.")
            return

        updated = 0
        for ds in datasets:
            qc = compute_quality_score(ds)
            ds.quality_score = qc["score"]
            ds.quality_status = qc["status"]
            ds.quality_issues = qc["issues"]
            updated += 1

        db.commit()
        print(f"Updated quality metrics for {updated} datasets.")
    finally:
        try:
            db.close()
        finally:
            try:
                next(db_gen, None)
            except Exception:
                pass


if __name__ == "__main__":
    main()

