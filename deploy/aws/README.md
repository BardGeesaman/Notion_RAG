# AWS Deployment (Lightsail + RDS) — Step-by-Step Guide

This guide deploys the Amprenta platform on AWS using:

- **Amazon RDS (PostgreSQL)** as the managed database
- **Amazon Lightsail** as a simple VM host for the application

It is intentionally **beginner-friendly** and “Phase-based”.

---

## Prerequisites

### Accounts & tools

- An **AWS account** with billing enabled
- **AWS CLI** installed locally
  - See AWS docs: [Installing the AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
- An SSH client (macOS/Linux `ssh`, Windows PowerShell/Windows Terminal)

### Credentials & access

- You should have **MFA enabled** and an **IAM user** (not using root access keys)
- You need these values ready for Phase 4:
  - `OPENAI_API_KEY`
  - `PINECONE_API_KEY`
  - Any other platform secrets used in `.env`

---

## Phase 1 — AWS Account Setup (MFA + IAM)

### 1. Enable MFA on the root account

1. AWS Console → **Account** / **Security credentials**
2. Enable MFA for **Root user**

**Screenshot reference**: “AWS Console → Security credentials → Multi-factor authentication (MFA)”.

### 2. Create an IAM admin user for day-to-day work

1. AWS Console → **IAM** → **Users** → “Create user”
2. Add the user to an **Admin** group (or attach `AdministratorAccess` policy)
3. Create **access keys** for CLI use (or use AWS SSO if your org prefers it)

**Screenshot reference**: “IAM → Users → Create user → Permissions”.

### 3. Configure AWS CLI locally

Run:

```bash
aws configure
```

Provide:
- **AWS Access Key ID**
- **AWS Secret Access Key**
- **Default region** (e.g. `us-east-1`)
- **Default output format** (e.g. `json`)

Verify:

```bash
aws sts get-caller-identity
```

---

## Phase 2 — RDS PostgreSQL Setup

### 1. Create the RDS instance

AWS Console → **RDS** → “Create database”:

- **Engine**: PostgreSQL
- **DB instance identifier**: `amprenta-prod`
- **Master username**: e.g. `amprenta`
- **Password**: generate and store securely
- **Public access**: **No** (recommended)
- **VPC**: default or a dedicated VPC
- **Backups**: enable (recommended)

**Screenshot reference**: “RDS → Create database → Connectivity”.

### 2. Create a security group rule allowing Lightsail to connect

You will later assign a **static IP** to your Lightsail instance.

In your RDS security group, add an inbound rule:
- **Type**: PostgreSQL
- **Port**: 5432
- **Source**: your Lightsail static IP (recommended) or the Lightsail instance security group (advanced)

### 3. Record connection details

You will need:
- RDS endpoint hostname (example: `amprenta-prod.xxxxx.us-east-1.rds.amazonaws.com`)
- Port (usually `5432`)
- Database name (e.g. `amprenta`)
- Username/password

### 4. (Optional) Connectivity test

From your local machine (if you have network path), test with:

```bash
psql "host=<RDS_HOST> port=5432 dbname=<DB_NAME> user=<DB_USER> sslmode=require"
```

If RDS is private-only, this test must be run from inside the VPC (e.g. from Lightsail or a bastion).

---

## Phase 3 — Lightsail Instance Setup

### 1. Create a Lightsail instance (Ubuntu)

AWS Console → **Lightsail** → “Create instance”:

- Platform: **Linux/Unix**
- Blueprint: **Ubuntu 22.04 LTS** (or 24.04 if preferred)
- Plan: choose based on expected usage

### 2. Attach a static IP (recommended)

Lightsail → **Networking** → “Create static IP” → attach to the instance.

You will use this IP for:
- SSH access
- RDS security group allowlist
- Optional DNS setup

### 3. Open required firewall ports

Lightsail “Networking” tab:
- **SSH**: 22 (restricted to your IP if possible)
- **HTTP**: 80 (if using a reverse proxy)
- **HTTPS**: 443 (if using TLS)
- **Streamlit**: 8501 (only if exposing directly; prefer proxy)
- **FastAPI**: 8000 (only if exposing directly; prefer proxy)

**Recommendation**: expose only 80/443 publicly and proxy internally.

### 4. SSH into the instance

Use the Lightsail browser SSH or local SSH:

```bash
ssh ubuntu@<LIGHTSAIL_STATIC_IP>
```

---

## Phase 4 — Application Deployment

### 1. Run the setup script

On the instance:

```bash
cd ~
git clone <REPO_URL> amprenta
cd amprenta

bash deploy/aws/setup-lightsail.sh
```

### 2. Create environment configuration (`.env`)

Create a `.env` in the repo root (example keys — adapt to your codebase):

```env
# Postgres (RDS)
POSTGRES_HOST=<RDS_HOST>
POSTGRES_PORT=5432
POSTGRES_DB=<DB_NAME>
POSTGRES_USER=<DB_USER>
POSTGRES_PASSWORD=<DB_PASSWORD>

# RAG / LLM
OPENAI_API_KEY=...
PINECONE_API_KEY=...

# Optional
DISABLE_AUTH=false
```

### 3. Run migrations

```bash
source ~/.bashrc
conda activate amprenta
alembic upgrade head
```

### 4. Start the app

#### Option A — Quick start (manual)

Streamlit:

```bash
conda activate amprenta
streamlit run streamlit_app/app.py --server.port 8501 --server.address 0.0.0.0
```

FastAPI (if needed):

```bash
conda activate amprenta
uvicorn amprenta_rag.api.main:app --host 0.0.0.0 --port 8000
```

#### Option B — Production start (systemd) (recommended)

Create systemd services to run Streamlit/FastAPI on boot. (This repo may already contain service templates; if not, create them in `/etc/systemd/system/`.)

**Screenshot reference**: “Lightsail → Instance → Networking → public URLs”.

---

## Troubleshooting

### “Cannot connect to RDS”

- Confirm Lightsail static IP is allowed in the **RDS security group** (port 5432).
- Confirm RDS is reachable from the Lightsail VPC/network path.
- Confirm `POSTGRES_HOST/USER/PASSWORD/DB` are correct.

### “ModuleNotFoundError” / dependency issues

- Ensure you activated the conda environment:
  - `conda activate amprenta`
- Re-run:
  - `pip install -r requirements.txt`

### “Streamlit is running but I can’t access it”

- Confirm Lightsail firewall includes port **8501** (or you are using a reverse proxy on 80/443).
- Confirm Streamlit is bound to `0.0.0.0`.

### “Pinecone/OpenAI auth errors”

- Verify `OPENAI_API_KEY` / `PINECONE_API_KEY` are set in `.env`.
- Restart the processes after updating `.env`.


