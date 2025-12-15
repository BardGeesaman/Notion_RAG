## Debug Protocol: Avoid Combinatorial Failures (One Variable at a Time)

When a feature “breaks,” it’s often multiple independent issues stacked together.
This protocol prevents combinatorial debugging by forcing isolation + deterministic smoke tests.

### Core rule
- **Change only one axis at a time**, then run the **smallest deterministic smoke test**.

### Axes (common independent variables)
- **Backend**: API routes wired? DB schema/data present? seed scripts?
- **Frontend**: browser console errors? widget manager assets? auth/cookies?
- **Notebook**: cell execution error? environment mismatch? missing package?
- **Container networking**: correct hostname (`host.docker.internal` vs `localhost` vs service name)?
- **Versions**: Voila/JupyterLab/ipywidgets/nbclient drift; stale assets.

### Deterministic smoke tests (examples)

#### API (SAR)
- **Goal**: prove data + computation work without UI.

```bash
curl -s http://localhost:8000/api/v1/sar/targets
curl -s http://localhost:8000/api/v1/sar/targets/CDK2/compounds | head
curl -s "http://localhost:8000/api/v1/sar/targets/CDK2/cliffs?similarity_threshold=0.6&fold_change=10&limit=1"
```

#### Voila UI (widgets rendering)
- **Goal**: prove widgets initialize, not just HTML loads.
- **Check**:
  - Page loads
  - widget controls appear (Text/Dropdown/Sliders)
  - browser console has no fatal module/version errors

### Stop-the-line heuristics
- If you don’t know whether the problem is **data** or **rendering**:
  - **Prove API first** (it’s faster, deterministic).
  - Then prove UI.
- If output is empty (`[]`) but “shouldn’t be”:
  - Verify the algorithm’s prerequisites (e.g., **RDKit can parse SMILES** → fingerprints exist) before tuning thresholds.

### “Last known good” discipline
- Once you find a working version combo, **pin it** (Dockerfile / requirements).
- Avoid “fix-forward by upgrading everything” unless you’re prepared to chase frontend asset mismatches.


