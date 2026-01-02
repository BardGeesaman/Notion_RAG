# P1 Auth Endpoint Remediation Plan

## Scope
- **230 endpoints** flagged as missing authentication
- **121 HIGH risk** (POST/PUT/DELETE/PATCH)
- **109 MEDIUM risk** (GET)

## Approach

For each endpoint:
1. **Review intent** - Is this intentionally public?
2. **Add auth dependency** - `current_user: User = Depends(get_current_user)`
3. **Add ownership check** - Use `can_view_entity()` / `can_edit_entity()` where applicable
4. **Test** - Verify endpoint requires auth

### Auth Pattern
```python
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.database.models import User
from fastapi import Depends

@router.get("/resource/{id}")
def get_resource(
    id: UUID,
    current_user: User = Depends(get_current_user),  # ADD THIS
    db: Session = Depends(get_db)
):
    # Optional: Add ownership check for user-owned resources
    resource = db.query(Resource).filter(Resource.id == id).first()
    if resource and resource.created_by != current_user.id:
        raise HTTPException(403, "Not authorized")
    return resource
```

### Intentionally Public Endpoints (skip auth)
- `/share_links/{token}/validate` - Public token validation
- `/entity_reviews/status-info` - Public status enum lookup
- Reference data lookups (e.g., `/phenotypes/{hpo_id}/genes`)

---

## Batch A: Core CRUD (22 endpoints) - Priority HIGH

### datasets.py (11 endpoints)
| Method | Path | Risk | Action |
|--------|------|------|--------|
| POST | / | HIGH | Add auth + created_by |
| GET | / | MED | Add auth + filter by user |
| GET | /{dataset_id} | MED | Add auth + ownership check |
| POST | /{dataset_id}/annotations | HIGH | Add auth + ownership |
| GET | /{dataset_id}/notebook | MED | Add auth + ownership |
| PATCH | /{dataset_id} | HIGH | Add auth + ownership |
| DELETE | /{dataset_id} | HIGH | Add auth + ownership |
| POST | /find | HIGH | Add auth |
| POST | /{dataset_id}/enrich | HIGH | Add auth + ownership |
| GET | /{dataset_id}/enrichment-status | MED | Add auth + ownership |
| POST | /{dataset_id}/expression-overlay | HIGH | Add auth + ownership |

### experiments.py (6 endpoints)
| Method | Path | Risk | Action |
|--------|------|------|--------|
| POST | / | HIGH | Add auth + created_by |
| GET | / | MED | Add auth + filter by user |
| GET | /{experiment_id} | MED | Add auth + ownership |
| POST | /{experiment_id}/annotations | HIGH | Add auth + ownership |
| PATCH | /{experiment_id} | HIGH | Add auth + ownership |
| DELETE | /{experiment_id} | HIGH | Add auth + ownership |

### programs.py (5 endpoints)
| Method | Path | Risk | Action |
|--------|------|------|--------|
| POST | / | HIGH | Add auth + created_by |
| GET | / | MED | Add auth + filter by user |
| GET | /{program_id} | MED | Add auth + ownership |
| PATCH | /{program_id} | HIGH | Add auth + ownership |
| DELETE | /{program_id} | HIGH | Add auth + ownership |

---

## Batch B: Research Data (29 endpoints) - Priority HIGH

### papers.py (10 endpoints)
### ip.py (9 endpoints)
### biophysical.py (10 endpoints)

---

## Batch C: Analysis (26 endpoints) - Priority MEDIUM

### planner.py (8 endpoints)
### single_cell.py (8 endpoints)
### crispr.py (6 endpoints)
### variants.py (6 endpoints)

---

## Batch D: ML/Compute (25 endpoints) - Priority MEDIUM

### notebook.py (6 endpoints)
### ml.py (4 endpoints)
### automl.py (4 endpoints)
### predictors.py (2 endpoints)
### scoring.py (3 endpoints)
### bayesian.py (3 endpoints)
### admet.py (2 endpoints)
### chemistry.py (2 endpoints)

---

## Batch E: Integrations (23 endpoints) - Priority MEDIUM

### connectivity.py (6 endpoints)
### graph.py (6 endpoints)
### subscriptions.py (6 endpoints)
### sync.py (5 endpoints)

---

## Batch F: Viz/Export (21 endpoints) - Priority LOW

### structures.py (6 endpoints)
### pathway_maps.py (5 endpoints)
### spectral.py (5 endpoints)
### ranking.py (5 endpoints)

---

## Batch G: Specialized (17 endpoints) - Priority LOW

### multi_omics.py (5 endpoints)
### multi_omics_viz.py (3 endpoints)
### poses.py (3 endpoints)
### docking.py (5 endpoints)
### biomarker.py (2 endpoints)

---

## Batch H: Remaining (~67 endpoints) - Priority LOW

### Files with <5 endpoints each:
- catalog.py (2)
- alerts.py (3)
- audit.py (4)
- jobs.py (3)
- export.py (4)
- features.py (5)
- projector.py (3)
- signatures.py (3)
- inline_annotations.py (1)
- compound_target.py (2)
- notebooks.py (4)
- viz3d.py (3)
- portfolio.py (4)
- explorer.py (4)
- batch.py (1)
- sphingolipid.py (1)
- qsar.py (3)
- sar.py (4)
- phenotypes.py (2)
- pockets.py (2)
- compounds.py (1)
- share_links.py (1 - SKIP, intentionally public)
- entity_reviews.py (1 - SKIP, status-info is public)

---

## Timeline
- Batch A: Day 1 (Core CRUD - highest impact)
- Batch B: Day 2-3 (Research data)
- Batch C: Day 3-4 (Analysis)
- Batch D-H: Day 5-10 (Lower priority)

## Success Criteria
- All 230 endpoints reviewed
- ~220+ endpoints secured with auth
- ~10 endpoints documented as intentionally public
- All existing tests pass
- New auth tests added for high-risk endpoints

