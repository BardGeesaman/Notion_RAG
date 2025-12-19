"""Seed synthetic HTS data (campaigns, plates, dose-response, hits)."""

from __future__ import annotations

import argparse
import random
from datetime import datetime, timedelta
from typing import Dict, List, Tuple

import numpy as np

from amprenta_rag.database.session import db_session
from amprenta_rag.models.chemistry import HTSCampaign, HTSResult, Compound


SIZE_PRESETS: Dict[str, Tuple[int, int]] = {
    "small": (1, 2),    # campaigns, plates
    "medium": (3, 10),
    "large": (10, 50),
}

WELLS_96 = [f"{chr(ord('A') + r)}{c:02d}" for r in range(8) for c in range(1, 13)]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed HTS demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo HTS data first.")
    parser.add_argument("--seed", type=int, default=2025)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> Tuple[int, int, int]:
    with db_session() as db:
        res_deleted = db.query(HTSResult).filter(HTSResult.result_id.like("DEMO_HTSRES_%")).delete(synchronize_session=False)
        camp_deleted = db.query(HTSCampaign).filter(HTSCampaign.campaign_id.like("DEMO_HTS_%")).delete(synchronize_session=False)
        cmp_deleted = db.query(Compound).filter(Compound.compound_id.like("DEMO_HTS_CMP_%")).delete(synchronize_session=False)
        db.commit()
    return camp_deleted, res_deleted, cmp_deleted


def _ensure_compounds(db, count: int, rng: random.Random) -> List[Compound]:
    needed = count
    created: List[Compound] = []
    existing = db.query(Compound).filter(Compound.compound_id.like("DEMO_HTS_CMP_%")).all()
    created.extend(existing)
    next_idx = len(existing) + 1
    while len(created) < needed:
        cid = f"DEMO_HTS_CMP_{next_idx:04d}"
        # Compound.smiles is non-nullable; use a tiny valid placeholder SMILES.
        c = Compound(compound_id=cid, smiles="C")
        db.add(c)
        created.append(c)
        next_idx += 1
    # Ensure PK defaults are populated (uuid) for any newly-added rows.
    db.flush()
    rng.shuffle(created)
    return created[:needed]


def _simulate_plate(compounds: List[Compound], plate_idx: int, rng: random.Random) -> List[Dict]:
    # Controls
    pos_wells = set(WELLS_96[:8])   # first 8 wells as positive
    neg_wells = set(WELLS_96[8:16]) # next 8 as negative

    results = []

    # precompute potencies for compounds (IC50 in uM)
    ic50_map = {c.compound_id: 10 ** rng.uniform(-2, 1) for c in compounds}  # between 0.01 and 10 uM

    # assign compounds to remaining wells
    sample_wells = [w for w in WELLS_96 if w not in pos_wells and w not in neg_wells]
    rng.shuffle(sample_wells)
    assignments = []
    for i, well in enumerate(sample_wells):
        compound = compounds[i % len(compounds)]
        # concentration range 0.01 - 30 uM
        conc = 10 ** rng.uniform(-2, 1.5)
        ic50 = ic50_map[compound.compound_id]
        hill = 1.0
        bottom = 0.05
        top = 1.0
        resp = bottom + (top - bottom) / (1 + (conc / ic50) ** hill)
        resp += rng.gauss(0, 0.05)
        resp = max(0.0, min(1.2, resp))
        assignments.append((well, compound, resp, conc, ic50))

    # compute plate stats for z-score
    raw_values = []
    for well in WELLS_96:
        if well in pos_wells:
            raw_values.append(0.05 + rng.gauss(0, 0.01))
        elif well in neg_wells:
            raw_values.append(1.0 + rng.gauss(0, 0.02))
        else:
            # placeholder; fill later
            raw_values.append(None)
    # replace placeholders with sample responses
    for idx, well in enumerate(WELLS_96):
        if raw_values[idx] is None:
            for a in assignments:
                if a[0] == well:
                    raw_values[idx] = a[2]
                    break

    mean = float(np.mean(raw_values))
    std = float(np.std(raw_values)) or 1.0

    # build results
    for well in WELLS_96:
        if well in pos_wells:
            raw = 0.05 + rng.gauss(0, 0.01)
            norm = raw
            z = (raw - mean) / std
            hit = True
            hit_cat = "control_positive"
            compound = None
        elif well in neg_wells:
            raw = 1.0 + rng.gauss(0, 0.02)
            norm = raw
            z = (raw - mean) / std
            hit = False
            hit_cat = "control_negative"
            compound = None
        else:
            assignment = next(a for a in assignments if a[0] == well)
            _, compound, raw, conc, ic50 = assignment
            norm = raw
            z = (raw - mean) / std
            hit = norm <= 0.3
            hit_cat = "hit" if hit else "non-hit"
        result = {
            "well": well,
            "raw": raw,
            "norm": norm,
            "z": z,
            "hit": hit,
            "hit_cat": hit_cat,
            "compound": compound,
        }
        results.append(result)
    return results


def _seed_hts(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int]:
    campaigns_count, plates_total = SIZE_PRESETS[size]
    rng = random.Random(seed_value)

    campaigns_created = results_created = 0

    with db_session() as db:
        # create campaigns
        campaigns: List[HTSCampaign] = []
        for idx in range(1, campaigns_count + 1):
            cid = f"DEMO_HTS_{idx:03d}"
            camp = db.query(HTSCampaign).filter(HTSCampaign.campaign_id == cid).first()
            if not camp:
                camp = HTSCampaign(
                    campaign_id=cid,
                    campaign_name=f"Demo HTS Campaign {idx}",
                    description="Synthetic HTS campaign",
                    assay_type="biochemical",
                    target="DEMO",
                    total_wells=96,
                    run_date=datetime.utcnow() - timedelta(days=idx),
                )
                db.add(camp)
                campaigns_created += 1
            campaigns.append(camp)

        # compounds pool
        compounds = _ensure_compounds(db, 200, rng)

        if not dry_run:
            db.commit()
            for c in campaigns:
                db.refresh(c)
        else:
            db.rollback()

    # assign plates across campaigns
    plates_per_camp = max(1, plates_total // campaigns_count)
    plate_counter = 1

    with db_session() as db:
        for camp in campaigns:
            for _ in range(plates_per_camp):
                if plate_counter > plates_total:
                    break
                # simulate plate
                plate_compounds = rng.sample(compounds, k=min(len(compounds), 60))
                plate_results = _simulate_plate(plate_compounds, plate_counter, rng)
                # create HTSResults
                for res_idx, res in enumerate(plate_results, start=1):
                    rid = f"DEMO_HTSRES_{plate_counter:03d}_{res_idx:03d}"
                    existing = db.query(HTSResult).filter(HTSResult.result_id == rid).first()
                    if existing:
                        continue
                    htsr = HTSResult(
                        result_id=rid,
                        campaign_id=camp.id,
                        compound_id=res["compound"].id if res["compound"] else compounds[0].id,
                        well_position=res["well"],
                        raw_value=res["raw"],
                        normalized_value=res["norm"],
                        z_score=res["z"],
                        hit_flag=res["hit"],
                        hit_category=res["hit_cat"],
                    )
                    db.add(htsr)
                    results_created += 1

                plate_counter += 1

            # update hit_count
            camp_results = db.query(HTSResult).filter(HTSResult.campaign_id == camp.id).all()
            camp.hit_count = sum(1 for r in camp_results if r.hit_flag)

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return campaigns_created, results_created


def main() -> None:
    args = _parse_args()
    if args.reset:
        camps, res, cmp = _reset_demo()
        print(f"Reset: campaigns={camps}, results={res}, compounds={cmp}")

    camps_created, res_created = _seed_hts(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Campaigns created: {camps_created}")
    print(f"Results created: {res_created}")


if __name__ == "__main__":
    main()

