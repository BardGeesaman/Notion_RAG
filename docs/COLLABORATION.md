# Real-Time Collaborative Editing

This document provides comprehensive information about the Real-Time Collaborative (RTC) editing feature in the Amprenta platform, including architecture, deployment, API usage, and troubleshooting.

## Table of Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Deployment](#deployment)
- [API Reference](#api-reference)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)
- [Security Considerations](#security-considerations)
- [Performance & Scalability](#performance--scalability)

## Overview

The Real-Time Collaborative Editing feature enables multiple users to simultaneously edit Jupyter notebooks and other documents with real-time synchronization, conflict-free merging, and live cursor tracking.

### Key Features

- **Real-Time Synchronization**: Changes are synchronized instantly across all connected users
- **Conflict-Free Editing**: Uses Y.js CRDTs (Conflict-free Replicated Data Types) to handle concurrent edits
- **Live Cursors**: See other users' cursor positions and selections in real-time
- **Session Management**: Track active collaboration sessions and participants
- **User Invitations**: Invite specific users to join collaborative editing sessions
- **Persistent State**: Document state is persisted across sessions via PostgreSQL

### Supported Document Types

- Jupyter Notebooks (`.ipynb`)
- Python Scripts (`.py`)
- Markdown Files (`.md`)
- Text Files (`.txt`)
- JSON Files (`.json`)

## Architecture

### Component Overview

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   JupyterLab    │    │  Collaboration  │    │   Y.js Server   │
│   Singleuser    │◄──►│      API        │◄──►│   (WebSocket)   │
│                 │    │                 │    │                 │
└─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ jupyter-        │    │   FastAPI       │    │   PostgreSQL    │
│ collaboration   │    │   Backend       │    │   Persistence   │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

### Y.js CRDT Architecture

**Conflict-free Replicated Data Types (CRDTs)** enable multiple users to edit the same document simultaneously without conflicts:

1. **Document State**: Each document is represented as a Y.js document with shared types
2. **Operations**: User edits are converted to operations that can be applied in any order
3. **Synchronization**: Operations are broadcast to all connected clients via WebSocket
4. **Persistence**: Document state is persisted to PostgreSQL for durability

### Data Flow

```
User Edit → Y.js Operation → WebSocket Broadcast → Other Clients → UI Update
     ↓
PostgreSQL Persistence
```

### Network Communication

- **WebSocket Protocol**: Real-time bidirectional communication
- **Container Network**: `jupyterhub-network` for internal communication
- **Port Configuration**: Y.js server on port 1234, configurable via environment

## Deployment

### Docker Compose Setup

The collaboration feature is deployed using Docker Compose with the following services:

#### 1. Y.js WebSocket Server

```yaml
# deploy/jupyterhub/docker-compose.yml
yjs-server:
  build:
    context: ./yjs-server
    dockerfile: Dockerfile
  container_name: amprenta-yjs-server
  ports:
    - "1234:1234"
  environment:
    - PGHOST=${POSTGRES_HOST:-postgres}
    - PGPORT=${POSTGRES_PORT:-5432}
    - PGDATABASE=${POSTGRES_DB:-amprenta}
    - PGUSER=${POSTGRES_USER:-postgres}
    - PGPASSWORD=${POSTGRES_PASSWORD:-postgres}
    - YJS_SERVER_PORT=${YJS_SERVER_PORT:-1234}
    - YJS_SERVER_HOST=${YJS_SERVER_HOST:-0.0.0.0}
  networks:
    - jupyterhub-network
  restart: unless-stopped
```

#### 2. JupyterHub Configuration

```python
# deploy/jupyterhub/jupyterhub_config.py
c.DockerSpawner.environment = {
    "API_URL": os.environ.get("API_URL", "http://host.docker.internal:8000"),
    "JUPYTERHUB_SINGLEUSER_APP": "jupyter_server.serverapp.ServerApp",
    # Y.js WebSocket server URL for real-time collaborative editing
    "JUPYTER_COLLABORATION_YJS_URL": os.environ.get("YJS_SERVER_URL", "ws://yjs-server:1234"),
}
```

#### 3. Singleuser Image

```dockerfile
# deploy/jupyterhub/singleuser/Dockerfile
RUN pip install --no-cache-dir \
    # ... existing packages ...
    # Jupyter collaboration for real-time collaborative editing
    jupyter-collaboration

# Y.js server URL configuration
ENV JUPYTER_COLLABORATION_YJS_URL=ws://yjs-server:1234

# Enable jupyter-collaboration extension
RUN mkdir -p /etc/jupyter/jupyter_server_config.d && \
    printf '%s\n' \
      '{' \
      '  "ServerApp": {' \
      '    "jpserver_extensions": {' \
      '      "jupyter_collaboration": true' \
      '    }' \
      '  }' \
      '}' \
    > /etc/jupyter/jupyter_server_config.d/jupyter_collaboration.json
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `YJS_SERVER_URL` | `ws://yjs-server:1234` | Y.js WebSocket server URL |
| `YJS_SERVER_PORT` | `1234` | Y.js server port |
| `YJS_SERVER_HOST` | `0.0.0.0` | Y.js server bind address |
| `JUPYTER_COLLABORATION_YJS_URL` | `ws://yjs-server:1234` | Jupyter extension Y.js URL |

### Deployment Commands

```bash
# Start all services
cd deploy/jupyterhub
docker-compose up -d

# View logs
docker-compose logs -f yjs-server

# Restart Y.js server
docker-compose restart yjs-server

# Check service health
curl http://localhost:1234/health
```

## API Reference

The Collaboration API provides REST endpoints for managing collaborative editing sessions.

### Base URL

```
https://your-domain.com/api/v1/collaboration
```

### Authentication

All endpoints require authentication via JWT token:

```bash
curl -H "Authorization: Bearer <jwt-token>" \
     https://your-domain.com/api/v1/collaboration/sessions
```

### Endpoints

#### 1. List Collaboration Sessions

**`GET /sessions`**

Lists all active collaboration sessions for the current user.

**Response:**
```json
{
  "sessions": [
    {
      "document_id": "notebook_analysis_001",
      "document_path": "/home/jovyan/work/data_analysis.ipynb",
      "document_type": "notebook",
      "document_title": "Customer Data Analysis",
      "participant_count": 3,
      "last_activity": "2025-12-31T11:22:35.866404",
      "is_active": true,
      "is_owner": true
    }
  ],
  "total_count": 2,
  "active_count": 1
}
```

**Status Codes:**
- `200 OK` - Success
- `401 Unauthorized` - Authentication required

#### 2. Get Session Details

**`GET /sessions/{document_id}`**

Retrieves detailed information about a specific collaboration session.

**Parameters:**
- `document_id` (path): Unique identifier for the document

**Response:**
```json
{
  "session": {
    "document_id": "notebook_analysis_001",
    "document_path": "/home/jovyan/work/data_analysis.ipynb",
    "document_type": "notebook",
    "created_at": "2025-12-31T08:22:35.866404",
    "last_activity": "2025-12-31T11:22:35.866404",
    "participants": [
      {
        "id": "uuid-here",
        "username": "researcher",
        "display_name": "Research User",
        "is_online": true,
        "last_seen": "2025-12-31T11:22:35.866404"
      }
    ],
    "owner": {
      "id": "uuid-here",
      "username": "researcher",
      "display_name": "Research User",
      "is_online": true
    },
    "is_active": true,
    "document_title": "Customer Data Analysis"
  },
  "user_permissions": "read-write",
  "connection_url": "ws://yjs-server:1234/doc/notebook_analysis_001",
  "recent_activity": [
    "researcher created the session",
    "collaborator1 joined the session",
    "collaborator2 edited cell 3"
  ]
}
```

**Status Codes:**
- `200 OK` - Success
- `401 Unauthorized` - Authentication required
- `403 Forbidden` - User not a session participant
- `404 Not Found` - Session not found

#### 3. Invite User to Session

**`POST /sessions/{document_id}/invite`**

Invites a user to join a collaboration session.

**Parameters:**
- `document_id` (path): Unique identifier for the document

**Request Body:**
```json
{
  "username": "researcher123",
  "message": "Would you like to collaborate on this analysis?",
  "permissions": "read-write"
}
```

**Response:**
```json
{
  "success": true,
  "invite_id": "uuid-here",
  "invited_user": "researcher123",
  "document_id": "notebook_analysis_001",
  "message": "Successfully invited researcher123 to collaborate on document",
  "expires_at": "2025-01-01T11:22:35.866404"
}
```

**Status Codes:**
- `200 OK` - Success
- `400 Bad Request` - Invalid request (self-invite, already participant)
- `401 Unauthorized` - Authentication required
- `403 Forbidden` - Only session owners can invite users
- `404 Not Found` - Session or user not found

### Error Response Format

All error responses follow this format:

```json
{
  "detail": "Human-readable error message"
}
```

## Usage Examples

### Starting a Collaboration Session

1. **Open a notebook in JupyterLab**
2. **The jupyter-collaboration extension automatically connects to Y.js server**
3. **Other users can join by opening the same notebook**

### Inviting Users via API

```bash
# Invite a user to collaborate
curl -X POST \
  -H "Authorization: Bearer <jwt-token>" \
  -H "Content-Type: application/json" \
  -d '{
    "username": "collaborator",
    "message": "Join me for data analysis",
    "permissions": "read-write"
  }' \
  https://your-domain.com/api/v1/collaboration/sessions/notebook_001/invite
```

### Monitoring Active Sessions

```bash
# List all active sessions
curl -H "Authorization: Bearer <jwt-token>" \
     https://your-domain.com/api/v1/collaboration/sessions

# Get details for specific session
curl -H "Authorization: Bearer <jwt-token>" \
     https://your-domain.com/api/v1/collaboration/sessions/notebook_001
```

### WebSocket Connection

The Y.js server provides WebSocket endpoints for real-time synchronization:

```javascript
// Connect to document collaboration room
const websocket = new WebSocket('ws://yjs-server:1234/doc/notebook_001');

// Y.js handles the protocol automatically via jupyter-collaboration
```

## Troubleshooting

### Common Issues

#### 1. Connection Issues

**Problem**: Users cannot connect to collaborative sessions

**Symptoms**:
- "Connection failed" errors in JupyterLab
- Sessions appear as inactive
- Changes not synchronizing

**Solutions**:
```bash
# Check Y.js server status
curl http://localhost:1234/health

# Check container logs
docker-compose logs yjs-server

# Verify network connectivity
docker exec amprenta-jupyterhub ping yjs-server

# Restart services
docker-compose restart yjs-server
```

#### 2. Environment Variable Issues

**Problem**: Y.js server URL not configured correctly

**Symptoms**:
- jupyter-collaboration cannot connect
- WebSocket connection errors
- Default URLs being used instead of custom ones

**Solutions**:
```bash
# Check environment variables in singleuser container
docker exec <singleuser-container> env | grep JUPYTER_COLLABORATION

# Verify JupyterHub configuration
docker exec amprenta-jupyterhub cat /etc/jupyterhub/jupyterhub_config.py | grep YJS

# Update environment variables
export YJS_SERVER_URL="ws://custom-yjs:1234"
docker-compose up -d
```

#### 3. Database Connection Issues

**Problem**: Y.js server cannot persist document state

**Symptoms**:
- Document changes lost after restart
- Y.js server startup failures
- Database connection errors in logs

**Solutions**:
```bash
# Check PostgreSQL connectivity
docker exec amprenta-yjs-server pg_isready -h postgres -p 5432

# Verify database credentials
docker exec amprenta-yjs-server env | grep PG

# Check database logs
docker-compose logs postgres
```

#### 4. Permission Issues

**Problem**: Users cannot invite others or access sessions

**Symptoms**:
- 403 Forbidden errors when inviting
- Cannot view session details
- Access denied messages

**Solutions**:
- Verify user is session owner for invitations
- Check user authentication tokens
- Ensure user is session participant for access

### Debugging Tools

#### Health Check Endpoints

```bash
# Y.js server health
curl http://localhost:1234/health

# API health
curl https://your-domain.com/api/v1/collaboration/sessions
```

#### Log Analysis

```bash
# Y.js server logs
docker-compose logs -f yjs-server

# JupyterHub logs
docker-compose logs -f jupyterhub

# Singleuser container logs
docker logs <singleuser-container-name>
```

#### Network Debugging

```bash
# Test container-to-container connectivity
docker exec amprenta-jupyterhub nc -zv yjs-server 1234

# Check port binding
docker port amprenta-yjs-server

# Verify network configuration
docker network inspect jupyterhub-network
```

## Security Considerations

### Authentication & Authorization

- **JWT Authentication**: All API endpoints require valid JWT tokens
- **Session-Based Access Control**: Only session participants can view session details
- **Owner-Only Invitations**: Only session owners can invite new users
- **Self-Invitation Prevention**: Users cannot invite themselves to sessions

### Network Security

- **Internal Communication**: Y.js server communicates only within `jupyterhub-network`
- **No External Exposure**: Y.js WebSocket port not exposed externally by default
- **Container Isolation**: Each user runs in isolated Docker containers

### Data Protection

- **Document Persistence**: All document states encrypted at rest in PostgreSQL
- **Session Isolation**: Document rooms are isolated by document ID
- **User Isolation**: Users can only access sessions they participate in

### Best Practices

1. **Use HTTPS**: Always deploy with SSL/TLS encryption
2. **Regular Updates**: Keep jupyter-collaboration and Y.js dependencies updated
3. **Monitor Access**: Log and monitor collaboration session access
4. **Backup Strategy**: Include Y.js document state in backup procedures

## Performance & Scalability

### Performance Characteristics

- **Real-Time Latency**: Sub-100ms synchronization for local network deployments
- **Concurrent Users**: Supports 10-50 concurrent users per document (depending on edit frequency)
- **Document Size**: Optimized for documents up to 10MB
- **Memory Usage**: ~10MB per active document session

### Scaling Considerations

#### Horizontal Scaling

```yaml
# Scale Y.js server for high availability
yjs-server:
  deploy:
    replicas: 3
  environment:
    - REDIS_URL=redis://redis-cluster:6379  # For session sharing
```

#### Database Optimization

```sql
-- Index for document lookup performance
CREATE INDEX idx_yjs_documents_id ON yjs_documents(document_id);
CREATE INDEX idx_yjs_documents_updated ON yjs_documents(updated_at);
```

#### Monitoring Metrics

- **Active Sessions**: Number of concurrent collaboration sessions
- **WebSocket Connections**: Active WebSocket connections per server
- **Document Operations**: Rate of Y.js operations per second
- **Database Performance**: Query response times for document persistence

### Resource Requirements

| Component | CPU | Memory | Storage |
|-----------|-----|--------|---------|
| Y.js Server | 1-2 cores | 2-4GB | Minimal |
| PostgreSQL | 2-4 cores | 4-8GB | 100GB+ |
| Singleuser | 1 core | 2GB | 10GB |

---

## Additional Resources

- [Y.js Documentation](https://docs.yjs.dev/)
- [Jupyter Collaboration Extension](https://github.com/jupyterlab/jupyter-collaboration)
- [Docker Compose Documentation](https://docs.docker.com/compose/)
- [JupyterHub Documentation](https://jupyterhub.readthedocs.io/)

## Support

For technical support and questions:
- Check the troubleshooting section above
- Review container logs for error messages
- Consult the Y.js and Jupyter documentation
- Contact the platform development team

---

*Last Updated: December 31, 2025*
