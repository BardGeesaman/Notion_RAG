## Amprenta Platform Deployment Guide

This guide describes how to deploy the **FastAPI** backend and **Streamlit** dashboard
for the Amprenta Multi-Omics Platform in a production-like environment:

- FastAPI (uvicorn) application: `amprenta_rag.api.main:app`
- Streamlit dashboard entrypoint: `scripts/run_dashboard.py`
- Database: PostgreSQL (system of record)

The examples assume a Linux host using **systemd** and **Nginx** as a reverse proxy.
Adapt paths and users to your environment.

---

### 1. Environment Setup

- **Clone the repository** to e.g. `/opt/amprenta`:

```bash
sudo mkdir -p /opt/amprenta
sudo chown "$USER":"$USER" /opt/amprenta
cd /opt/amprenta
git clone <YOUR_REPO_URL> .
```

- **Create virtual environment** (recommended):

```bash
cd /opt/amprenta
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

- **Environment variables**: create `/opt/amprenta/.env` (or use systemd `EnvironmentFile`):

```bash
# /opt/amprenta/.env

# Core
AM_PRENTA_ENV=production
PYTHONUNBUFFERED=1

# FastAPI / UVicorn
API_HOST=0.0.0.0
API_PORT=8000
API_RELOAD=0

# Streamlit
DASHBOARD_HOST=0.0.0.0
DASHBOARD_PORT=8501

# Database (example)
POSTGRES_HOST=127.0.0.1
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=amprenta
POSTGRES_PASSWORD=change_me
```

Ensure the user running the services has read access to `.env`.

---

### 2. Database Initialization

Assuming PostgreSQL is installed and reachable:

1. **Create database and user** (example):

```bash
sudo -u postgres psql <<'SQL'
CREATE USER amprenta WITH PASSWORD 'change_me';
CREATE DATABASE amprenta OWNER amprenta;
GRANT ALL PRIVILEGES ON DATABASE amprenta TO amprenta;
SQL
```

2. **Run migrations** (from repo root):

```bash
cd /opt/amprenta
source .venv/bin/activate
alembic upgrade head
```

3. Optionally, load **seed data** as needed using provided scripts under `scripts/` or documented ingestion commands.

---

### 3. Local Process Management (Dev / Simple Deploy)

From the repository root:

```bash
cd /opt/amprenta
source .venv/bin/activate

# Start API (FastAPI + uvicorn)
bash scripts/start_api.sh

# In a separate shell, start dashboard
bash scripts/start_dashboard.sh
```

Logs:

```bash
tail -f logs/api.log
tail -f logs/dashboard.log
```

---

### 4. Systemd Service Files

The following examples assume:

- Repo root: `/opt/amprenta`
- Virtualenv: `/opt/amprenta/.venv`
- Service user: `amprenta`

Adjust paths/users as needed.

#### 4.1 FastAPI Service (`amprenta-api.service`)

Create `/etc/systemd/system/amprenta-api.service`:

```ini
[Unit]
Description=Amprenta FastAPI Service
After=network.target

[Service]
Type=simple
User=amprenta
WorkingDirectory=/opt/amprenta
EnvironmentFile=/opt/amprenta/.env
ExecStart=/bin/bash -lc 'source /opt/amprenta/.venv/bin/activate && bash scripts/start_api.sh'
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

#### 4.2 Streamlit Dashboard Service (`amprenta-dashboard.service`)

Create `/etc/systemd/system/amprenta-dashboard.service`:

```ini
[Unit]
Description=Amprenta Streamlit Dashboard
After=network.target

[Service]
Type=simple
User=amprenta
WorkingDirectory=/opt/amprenta
EnvironmentFile=/opt/amprenta/.env
ExecStart=/bin/bash -lc 'source /opt/amprenta/.venv/bin/activate && bash scripts/start_dashboard.sh'
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

#### 4.3 Enable and Manage Services

```bash
sudo systemctl daemon-reload
sudo systemctl enable amprenta-api.service amprenta-dashboard.service
sudo systemctl start amprenta-api.service amprenta-dashboard.service

# Check status
sudo systemctl status amprenta-api.service
sudo systemctl status amprenta-dashboard.service
```

Logs are written to:

- `logs/api.log`
- `logs/dashboard.log`

and also visible via `journalctl -u amprenta-api.service` etc.

---

### 5. Nginx Reverse Proxy Configuration

Example Nginx server block to expose:

- FastAPI under `/api/`
- Streamlit dashboard at `/`

Assume services listen on:

- FastAPI: `127.0.0.1:8000`
- Dashboard: `127.0.0.1:8501`

Create `/etc/nginx/sites-available/amprenta.conf`:

```nginx
server {
    listen 80;
    server_name amprenta.example.com;

    # Redirect HTTP to HTTPS (optional; if TLS is configured separately)
    # return 301 https://$host$request_uri;

    # Proxy to Streamlit dashboard (root)
    location / {
        proxy_pass http://127.0.0.1:8501/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Proxy to FastAPI under /api/
    location /api/ {
        proxy_pass http://127.0.0.1:8000/api/;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Optional: health check endpoint
    location /healthz {
        proxy_pass http://127.0.0.1:8000/health;
    }
}
```

Enable the site and reload Nginx:

```bash
sudo ln -s /etc/nginx/sites-available/amprenta.conf /etc/nginx/sites-enabled/amprenta.conf
sudo nginx -t
sudo systemctl reload nginx
```

For HTTPS, use **Certbot** or your TLS termination method of choice and update `listen`/`server_name` accordingly.

---

### 6. Operational Checklist

- **On first deploy**
  - [ ] Repo cloned to `/opt/amprenta`
  - [ ] Virtualenv created and dependencies installed
  - [ ] `.env` configured (API/Dashboard ports, DB credentials, secrets)
  - [ ] Postgres database created and migrations run (`alembic upgrade head`)
  - [ ] `logs/` directory exists and writable by service user
  - [ ] Systemd services enabled and running
  - [ ] Nginx config live and passing health checks

- **To update code**
  - [ ] `cd /opt/amprenta && git pull`
  - [ ] `source .venv/bin/activate && pip install -r requirements.txt` (if deps changed)
  - [ ] `alembic upgrade head` (if migrations added)
  - [ ] `sudo systemctl restart amprenta-api.service amprenta-dashboard.service`

---

### 7. Troubleshooting

- **API not responding**:
  - Check `logs/api.log`
  - `sudo systemctl status amprenta-api.service`
  - Confirm DB connectivity and configuration in `.env`

- **Dashboard not loading**:
  - Check `logs/dashboard.log`
  - `sudo systemctl status amprenta-dashboard.service`
  - Ensure Streamlit and all dashboard dependencies are installed in the venv

- **Nginx 502/504 errors**:
  - Ensure backend services are listening on expected ports (8000/8501)
  - Check Nginx error logs (e.g. `/var/log/nginx/error.log`)


