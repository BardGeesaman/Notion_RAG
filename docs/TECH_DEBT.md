# Technical Debt Tracker

Items identified during code reviews that are non-blocking but should be addressed.

## Priority Levels
- **P2**: Should fix before production
- **P3**: Nice to have / future enhancement

---

## Open Items

### Repository Health
| ID | Priority | Description | File | Added |
|----|----------|-------------|------|-------|
| HEALTH-1 | P3 | Add sync event tracking API for history chart | repository_health.py | 2024-12-18 |

---

## Resolved Items

| ID | Resolution | Commit | Date |
|----|------------|--------|------|
| API-1 | Fixed datasets 500 error (None→[] validator) | pending | 2024-12-18 |
| JOB-1 | Keyword filtering implemented | pending | 2024-12-18 |
| JOB-2 | Notification delivery (MVP logging) | pending | 2024-12-18 |
| CAT-3 | Format dates (e.g., "2 days ago") | pending | 2024-12-18 |
| NAV-2 | Auto-mark alert read on View navigation | pending | 2024-12-18 |
| BELL-1 | Show dataset title instead of UUID | pending | 2024-12-18 |
| DS-3 | Distinguish 404 vs 500 errors in UI | pending | 2024-12-18 |
| SUB-1 | Timestamps added to SubscriptionResponse | completed | 2025-12-19 |
| DS-2 | UUID validation implemented | completed | 2025-12-19 |
| BELL-2 | Confirmation dialog for "mark all read" | completed | 2025-12-19 |
| BELL-3 | Navigation link to dataset from alert | completed | 2025-12-19 |
| CAT-1 | Query optimization implemented | completed | 2025-12-19 |
| CAT-2 | Pagination buttons added | completed | 2025-12-19 |
| CAT-4 | Response schemas created | completed | 2025-12-19 |
| DS-1 | "Show more" pagination for features | completed | 2025-12-19 |
| NAV-1 | Breadcrumb trail navigation | completed | 2025-12-19 |
| HEALTH-2 | Configurable health score | completed | 2025-12-19 |
| HEALTH-3 | Auto-refresh timer | completed | 2025-12-19 |
| HEALTH-4 | Sync frequency and success rate metrics | completed | 2025-12-19 |

