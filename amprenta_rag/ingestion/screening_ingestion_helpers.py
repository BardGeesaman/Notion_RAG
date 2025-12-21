from __future__ import annotations

from typing import List

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Compound, HTSCampaign, HTSResult, BiochemicalResult


def insert_compound_pg(compound: Compound) -> None:
    with db_session() as db:
        obj = Compound(
            compound_id=compound.compound_id,
            smiles=compound.smiles,
            inchi_key=getattr(compound, "inchi_key", None),
            canonical_smiles=getattr(compound, "canonical_smiles", None),
            molecular_formula=getattr(compound, "molecular_formula", None),
            molecular_weight=getattr(compound, "molecular_weight", None),
            logp=getattr(compound, "logp", None),
            hbd_count=getattr(compound, "hbd_count", None),
            hba_count=getattr(compound, "hba_count", None),
            rotatable_bonds=getattr(compound, "rotatable_bonds", None),
        )
        db.add(obj)
        db.commit()


def insert_hts_campaign_pg(campaign: HTSCampaign) -> None:
    with db_session() as db:
        obj = HTSCampaign(
            campaign_id=campaign.campaign_id,
            campaign_name=campaign.campaign_name,
            description=getattr(campaign, "description", None),
            assay_type=getattr(campaign, "assay_type", None),
            target=getattr(campaign, "target", None),
            library_id=getattr(campaign, "library_id", None),
            total_wells=getattr(campaign, "total_wells", None),
            hit_count=getattr(campaign, "hit_count", None),
            run_date=getattr(campaign, "run_date", None),
        )
        db.add(obj)
        db.commit()


def insert_hts_results_pg(results: List[HTSResult]) -> None:
    if not results:
        return
    with db_session() as db:
        campaign_ids = {str(r.campaign_id) for r in results if getattr(r, "campaign_id", None)}
        campaign_map = {}
        if campaign_ids:
            db_campaigns = db.query(HTSCampaign).filter(HTSCampaign.campaign_id.in_(campaign_ids)).all()
            campaign_map = {str(c.campaign_id): c.id for c in db_campaigns}

        compound_ids = {str(r.compound_id) for r in results if getattr(r, "compound_id", None)}
        compound_map = {}
        if compound_ids:
            db_compounds = db.query(Compound).filter(Compound.compound_id.in_(compound_ids)).all()
            compound_map = {str(c.compound_id): c.id for c in db_compounds}

        objects = []
        for r in results:
            campaign_uuid = campaign_map.get(str(r.campaign_id)) if getattr(r, "campaign_id", None) else None
            compound_uuid = compound_map.get(str(r.compound_id)) if getattr(r, "compound_id", None) else None
            obj = HTSResult(
                result_id=r.result_id,
                campaign_id=campaign_uuid,
                compound_id=compound_uuid,
                well_position=getattr(r, "well_position", None),
                raw_value=getattr(r, "raw_value", None),
                normalized_value=getattr(r, "normalized_value", None),
                z_score=getattr(r, "z_score", None),
                hit_flag=getattr(r, "hit_flag", None),
                hit_category=getattr(r, "hit_category", None),
            )
            objects.append(obj)

        if objects:
            db.bulk_save_objects(objects)
            db.commit()


def insert_biochemical_results_pg(results: List[BiochemicalResult]) -> None:
    if not results:
        return
    with db_session() as db:
        compound_ids = {str(r.compound_id) for r in results if getattr(r, "compound_id", None)}
        compound_map = {}
        if compound_ids:
            db_compounds = db.query(Compound).filter(Compound.compound_id.in_(compound_ids)).all()
            compound_map = {str(c.compound_id): c.id for c in db_compounds}

        objects = []
        for r in results:
            compound_uuid = compound_map.get(str(r.compound_id)) if getattr(r, "compound_id", None) else None
            obj = BiochemicalResult(
                result_id=r.result_id,
                compound_id=compound_uuid,
                assay_name=getattr(r, "assay_name", None),
                target=getattr(r, "target", None),
                ic50=getattr(r, "ic50", None),
                ec50=getattr(r, "ec50", None),
                ki=getattr(r, "ki", None),
                kd=getattr(r, "kd", None),
                activity_type=getattr(r, "activity_type", None),
                units=getattr(r, "units", None),
                run_date=getattr(r, "run_date", None),
            )
            objects.append(obj)

        if objects:
            db.bulk_save_objects(objects)
            db.commit()

