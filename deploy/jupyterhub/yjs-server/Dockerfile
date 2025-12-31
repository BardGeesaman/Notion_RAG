# Y.js WebSocket Server Dockerfile
FROM node:18-alpine

# Set working directory
WORKDIR /app

# Install system dependencies for native modules
RUN apk add --no-cache \
    python3 \
    make \
    g++ \
    postgresql-client

# Copy package files
COPY package*.json ./

# Install Node.js dependencies
RUN npm ci --only=production && \
    npm cache clean --force

# Copy application code
COPY server.js ./

# Create non-root user
RUN addgroup -g 1001 -S nodejs && \
    adduser -S yjs -u 1001 -G nodejs

# Change ownership of app directory
RUN chown -R yjs:nodejs /app

# Switch to non-root user
USER yjs

# Health check
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD wget --no-verbose --tries=1 --spider http://localhost:1234/health || exit 1

# Expose port
EXPOSE 1234

# Environment variables
ENV NODE_ENV=production
ENV PORT=1234
ENV HOST=0.0.0.0

# Start command
CMD ["node", "server.js"]