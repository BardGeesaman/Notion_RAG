# Technical Debt Tracker

Items identified during code reviews that are non-blocking but should be addressed.

## Priority Levels
- **P2**: Should fix before production
- **P3**: Nice to have / future enhancement

---

## Open Items

### External Catalog (997e198)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| CAT-1 | P3 | Optimize double query (count + all) | catalog_service.py | 2024-12-18 |
| CAT-2 | P3 | Add prev/next pagination buttons | external_catalog.py | 2024-12-18 |
| CAT-3 | P3 | Format dates (e.g., "2 days ago") | external_catalog.py | 2024-12-18 |
| CAT-4 | P3 | Create proper response schemas | catalog.py | 2024-12-18 |

### Subscriptions (09cb143)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| SUB-1 | P3 | Add timestamps to SubscriptionResponse | schemas.py | 2024-12-18 |

### Background Job (8a78a11)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| — | — | — | — | — |

### Alerts Bell (5e5ddea)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| BELL-1 | P3 | Show dataset title instead of UUID | alerts_bell.py | 2024-12-18 |
| BELL-2 | P3 | Add confirmation for "mark all read" | alerts_bell.py | 2024-12-18 |
| BELL-3 | P3 | Add navigation link to dataset from alert | alerts_bell.py | 2024-12-18 |

### Dataset Details (pending)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| DS-1 | P3 | Add "Show more" pagination for features | dataset_details.py | 2024-12-18 |
| DS-2 | P3 | Validate UUID format client-side | dataset_details.py | 2024-12-18 |
| DS-3 | P3 | Distinguish 404 vs 500 errors in UI | dataset_details.py | 2024-12-18 |

### Repository Health (pending commit)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| HEALTH-1 | P3 | Add sync event tracking API for history chart | repository_health.py | 2024-12-18 |
| HEALTH-2 | P3 | Make health score configurable/formula-based | repository_health.py | 2024-12-18 |
| HEALTH-3 | P3 | Add auto-refresh timer for live updates | repository_health.py | 2024-12-18 |
| HEALTH-4 | P3 | Add sync frequency and success rate metrics | repository_health.py | 2024-12-18 |

### Navigation Links (pending commit)
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| NAV-1 | P3 | Add breadcrumb trail for navigation history | dataset_details.py | 2024-12-18 |
| NAV-2 | P3 | Auto-mark alert read on View navigation | alerts_bell.py | 2024-12-18 |

---

## Resolved Items

| ID | Resolution | Commit | Date |
|----|------------|--------|------|
| API-1 | Fixed datasets 500 error (None→[] validator) | pending | 2024-12-18 |
| JOB-1 | Keyword filtering implemented | pending | 2024-12-18 |
| JOB-2 | Notification delivery (MVP logging) | pending | 2024-12-18 |

