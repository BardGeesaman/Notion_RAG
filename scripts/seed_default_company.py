from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
from typing import Iterable, Optional

from amprenta_rag.database.models import Company, User
from amprenta_rag.database.session import db_session


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Seed a default Company and assign existing users to it.")
    ap.add_argument("--name", type=str, default="Default Organization")
    ap.add_argument("--subdomain", type=str, default="default")
    ap.add_argument("--trial-days", type=int, default=14)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(list(argv) if argv is not None else None)

    trial_ends = datetime.now(timezone.utc) + timedelta(days=int(args.trial_days))

    with db_session() as db:
        comp = db.query(Company).filter(Company.subdomain == args.subdomain).first()
        if comp is None:
            comp = Company(
                name=str(args.name),
                subdomain=str(args.subdomain),
                status="active",
                trial_ends_at=trial_ends,
            )
            db.add(comp)
            if not args.dry_run:
                db.commit()
                db.refresh(comp)

        users = db.query(User).all()
        updated = 0
        for u in users:
            if getattr(u, "company_id", None) is None:
                u.company_id = comp.id
                updated += 1
            if getattr(u, "company_role", None) is None:
                u.company_role = "member"
        if not args.dry_run:
            db.commit()

    print(f"Default company id={comp.id} subdomain={comp.subdomain} updated_users={updated} dry_run={args.dry_run}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


