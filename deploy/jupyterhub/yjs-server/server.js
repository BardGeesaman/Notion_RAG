#!/usr/bin/env node

/**
 * Y.js WebSocket Server with PostgreSQL Persistence
 * 
 * Provides collaborative document editing via Y.js CRDTs with persistent storage
 * in PostgreSQL for the Amprenta platform.
 */

require('dotenv').config();

const express = require('express');
const http = require('http');
const WebSocket = require('ws');
const Y = require('yjs');
const { setupWSConnection } = require('y-websocket/bin/utils');
const { Pool } = require('pg');

// Configuration
const PORT = process.env.PORT || 1234;
const HOST = process.env.HOST || '0.0.0.0';

// PostgreSQL connection
const pool = new Pool({
  host: process.env.POSTGRES_HOST || 'postgres',
  port: process.env.POSTGRES_PORT || 5432,
  database: process.env.POSTGRES_DB || 'amprenta',
  user: process.env.POSTGRES_USER || 'postgres',
  password: process.env.POSTGRES_PASSWORD || 'password',
  max: 20,
  idleTimeoutMillis: 30000,
  connectionTimeoutMillis: 2000,
});

// Document persistence
class PostgresPersistence {
  constructor(pool) {
    this.pool = pool;
    this.initTables();
  }

  async initTables() {
    try {
      await this.pool.query(`
        CREATE TABLE IF NOT EXISTS yjs_documents (
          name VARCHAR(255) PRIMARY KEY,
          state BYTEA NOT NULL,
          created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
          updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        )
      `);
      console.log('[YJS-SERVER] PostgreSQL tables initialized');
    } catch (error) {
      console.error('[YJS-SERVER] Failed to initialize tables:', error);
    }
  }

  async getYDoc(docname) {
    try {
      const result = await this.pool.query(
        'SELECT state FROM yjs_documents WHERE name = $1',
        [docname]
      );
      
      if (result.rows.length > 0) {
        const ydoc = new Y.Doc();
        const state = new Uint8Array(result.rows[0].state);
        Y.applyUpdate(ydoc, state);
        console.log(`[YJS-SERVER] Loaded document: ${docname}`);
        return ydoc;
      } else {
        console.log(`[YJS-SERVER] Creating new document: ${docname}`);
        return new Y.Doc();
      }
    } catch (error) {
      console.error(`[YJS-SERVER] Error loading document ${docname}:`, error);
      return new Y.Doc();
    }
  }

  async storeUpdate(docname, update) {
    try {
      await this.pool.query(`
        INSERT INTO yjs_documents (name, state, updated_at)
        VALUES ($1, $2, NOW())
        ON CONFLICT (name)
        DO UPDATE SET state = $2, updated_at = NOW()
      `, [docname, Buffer.from(update)]);
      
      console.log(`[YJS-SERVER] Stored update for document: ${docname}`);
    } catch (error) {
      console.error(`[YJS-SERVER] Error storing document ${docname}:`, error);
    }
  }

  async writeState(docname, ydoc) {
    const state = Y.encodeStateAsUpdate(ydoc);
    await this.storeUpdate(docname, state);
  }
}

// Initialize persistence
const persistence = new PostgresPersistence(pool);

// Document storage
const docs = new Map();

// Get or create document
const getYDoc = async (docname, gc = true) => {
  if (docs.has(docname)) {
    return docs.get(docname);
  }
  
  const doc = await persistence.getYDoc(docname);
  if (gc) {
    doc.gc = gc;
  }
  
  docs.set(docname, doc);
  
  // Set up persistence on document updates
  doc.on('update', async (update) => {
    await persistence.storeUpdate(docname, update);
  });
  
  return doc;
};

// Express app for health checks
const app = express();
app.use(express.json());

// Health check endpoint
app.get('/health', (req, res) => {
  res.json({ 
    status: 'healthy',
    service: 'yjs-websocket-server',
    timestamp: new Date().toISOString(),
    documents: docs.size
  });
});

// Metrics endpoint
app.get('/metrics', (req, res) => {
  res.json({
    documents_count: docs.size,
    active_documents: Array.from(docs.keys()),
    uptime: process.uptime(),
    memory_usage: process.memoryUsage()
  });
});

// Create HTTP server
const server = http.createServer(app);

// WebSocket server
const wss = new WebSocket.Server({ 
  server,
  path: '/yjs'
});

wss.on('connection', (conn, req) => {
  console.log('[YJS-SERVER] New WebSocket connection');
  
  // Set up Y.js connection with custom document getter
  setupWSConnection(conn, req, { 
    docName: req.url?.slice(1) || 'default',
    getYDoc: getYDoc
  });
});

// Graceful shutdown
const gracefulShutdown = async (signal) => {
  console.log(`[YJS-SERVER] Received ${signal}, shutting down gracefully...`);
  
  // Save all documents
  for (const [docname, doc] of docs.entries()) {
    try {
      await persistence.writeState(docname, doc);
      console.log(`[YJS-SERVER] Saved document: ${docname}`);
    } catch (error) {
      console.error(`[YJS-SERVER] Error saving document ${docname}:`, error);
    }
  }
  
  // Close connections
  wss.close(() => {
    console.log('[YJS-SERVER] WebSocket server closed');
    pool.end(() => {
      console.log('[YJS-SERVER] PostgreSQL pool closed');
      process.exit(0);
    });
  });
};

// Handle shutdown signals
process.on('SIGTERM', () => gracefulShutdown('SIGTERM'));
process.on('SIGINT', () => gracefulShutdown('SIGINT'));

// Start server
server.listen(PORT, HOST, () => {
  console.log(`[YJS-SERVER] Server running on ${HOST}:${PORT}`);
  console.log(`[YJS-SERVER] WebSocket endpoint: ws://${HOST}:${PORT}/yjs`);
  console.log(`[YJS-SERVER] Health check: http://${HOST}:${PORT}/health`);
});

// Handle uncaught exceptions
process.on('uncaughtException', (error) => {
  console.error('[YJS-SERVER] Uncaught exception:', error);
  gracefulShutdown('uncaughtException');
});

process.on('unhandledRejection', (reason, promise) => {
  console.error('[YJS-SERVER] Unhandled rejection at:', promise, 'reason:', reason);
  gracefulShutdown('unhandledRejection');
});