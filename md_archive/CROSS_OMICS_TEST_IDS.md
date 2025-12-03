# Cross-Omics Test IDs - Ready to Use!

## ✅ Available Test Data

### 🔖 Signatures

1. **ALS-CSF-Core-6Ceramides**
   - Short ID: `ALS-CSF-Core-6Cer-v1`
   - Page ID: `18eb23f2-ceec-45ed-a19e-9b540b85922d`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --signature 18eb23f2-ceec-45ed-a19e-9b540b85922d
     ```

2. **test_signature_verification**
   - Short ID: `test-signature-verification`
   - Page ID: `2beadf61-42ab-81cc-a119-f11b12a611bd`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --signature 2beadf61-42ab-81cc-a119-f11b12a611bd
     ```

### 🧬 Gene Features

1. **ALB**
   - Page ID: `2beadf61-42ab-8141-aff1-f0970b9424f1`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "gene:ALB"
     ```

2. **GFAP**
   - Page ID: `2beadf61-42ab-8190-b023-cb98ceb0fd99`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "gene:GFAP"
     ```

3. **TUBB3**
   - Page ID: `2beadf61-42ab-81b9-9387-c58253b52599`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "gene:TUBB3"
     ```

4. **APOE**
   - Page ID: `2beadf61-42ab-81d2-9750-cd40168c60ae`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "gene:APOE"
     ```

### 🧬 Protein Features

1. **GFAP**
   - Page ID: `2beadf61-42ab-8120-8c14-dd7481c93b09`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "protein:GFAP"
     ```

2. **APOE**
   - Page ID: `2beadf61-42ab-8138-921a-f1404956f85a`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "protein:APOE"
     ```

3. **TUBB3**
   - Page ID: `2beadf61-42ab-81cc-8c55-df77faefc28f`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "protein:TUBB3"
     ```

4. **ALB**
   - Page ID: `2beadf61-42ab-81e3-8adc-d720711cb36b`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "protein:ALB"
     ```

### 🧬 Metabolite Features

1. **Glutamine**
   - Page ID: `2beadf61-42ab-8152-b5a5-c9e9496bc347`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "metabolite:Glutamine"
     ```

2. **Glutamate** ✅ (Already tested - works!)
   - Page ID: `2beadf61-42ab-816f-bf28-c8c4ed8a62e7`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "metabolite:Glutamate"
     ```

3. **Serine**
   - Page ID: `2beadf61-42ab-8181-b167-fb2b1f9ad00e`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "metabolite:Serine"
     ```

4. **Lactate**
   - Page ID: `2beadf61-42ab-81a3-9d60-fc6f49f81e6c`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "metabolite:Lactate"
     ```

### 🧬 Lipid Features

1. **Cer(d18:1/18:0)**
   - Page ID: `2beadf61-42ab-8123-9a72-f791b485f632`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "lipid:Cer(d18:1/18:0)"
     ```

2. **SM(d18:1/24:1)**
   - Page ID: `2beadf61-42ab-812a-ae03-fa28c098438e`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "lipid:SM(d18:1/24:1)"
     ```

3. **Cer(d18:1/24:0)**
   - Page ID: `2beadf61-42ab-8163-9572-c1970a29cf24`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "lipid:Cer(d18:1/24:0)"
     ```

4. **Cer(d18:1/22:0)**
   - Page ID: `2beadf61-42ab-8163-a1ae-cd3fb1668263`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "lipid:Cer(d18:1/22:0)"
     ```

5. **HexCer(d18:1/22:0)**
   - Page ID: `2beadf61-42ab-8166-97f7-d1a104dc6301`
   - **Test Command:**
     ```bash
     python test_cross_omics.py --feature "lipid:HexCer(d18:1/22:0)"
     ```

---

## 🧪 Experiments

✅ **Experiments database is configured and working!**

1. **🧬 Multi-Omics Experiment**
   - Page ID: `288fa0da-a8a5-49ae-a5ee-2d5b94eaf661`
   - **Test Command:**
     ```bash
     # Note: Cross-omics experiment summary not yet implemented
     # But experiment data is available for cross-omics program summaries
     ```

2. **+6 more experiments** - Run `list_cross_omics_test_ids.py` to see all

## 📊 Programs & Datasets

⚠️ **Programs Database:** Currently has connection issues (HTTP 400). Needs troubleshooting.
✅ **Datasets Database:** Configured and working via `NOTION_EXP_DATA_DB_ID`

To list:
```bash
python scripts/list_cross_omics_test_ids.py --type experiments --limit 10
python scripts/list_cross_omics_test_ids.py --type datasets --limit 10
```

---

## 🚀 Quick Test Commands

### Test All Feature Types:

```bash
# Metabolites (already tested Glutamate ✅)
python test_cross_omics.py --feature "metabolite:Glutamine"
python test_cross_omics.py --feature "metabolite:Serine"
python test_cross_omics.py --feature "metabolite:Lactate"

# Genes
python test_cross_omics.py --feature "gene:GFAP"
python test_cross_omics.py --feature "gene:APOE"

# Proteins
python test_cross_omics.py --feature "protein:GFAP"
python test_cross_omics.py --feature "protein:APOE"

# Lipids
python test_cross_omics.py --feature "lipid:Cer(d18:1/18:0)"
python test_cross_omics.py --feature "lipid:SM(d18:1/24:1)"
```

### Test Signatures:

```bash
# ALS signature
python test_cross_omics.py --signature 18eb23f2-ceec-45ed-a19e-9b540b85922d

# Test signature
python test_cross_omics.py --signature 2beadf61-42ab-81cc-a119-f11b12a611bd
```

---

## 📝 Refresh Test IDs

To get the latest list of available IDs:

```bash
python scripts/list_cross_omics_test_ids.py --type all --limit 10
```

---

**Ready to test! Start with any feature or signature above!** 🎉

