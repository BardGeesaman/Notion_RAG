#!/bin/bash
# Setup Raycast Script Commands for Amprenta RAG

set -e

SCRIPT_DIR="$HOME/.raycast-rag"
RAG_ROOT="$HOME/Documents/Notion RAG"

echo "Creating Raycast script directory at: $SCRIPT_DIR"
mkdir -p "$SCRIPT_DIR"

#############################
# 1. Reindex Entire RAG
#############################
cat > "$SCRIPT_DIR/rag_reindex_all.sh" << 'EOF'
#!/bin/bash
# @raycast.title Reindex Entire RAG System
# @raycast.author Amprenta
# @raycast.authorURL https://amprenta.com
# @raycast.description Rebuild Zotero collection and ingest emails into the RAG system.
# @raycast.icon 🔄
# @raycast.mode fullOutput
# @raycast.package Amprenta RAG
# @raycast.refreshTime 1h

RAG_ROOT="$HOME/Documents/Notion RAG"
cd "$RAG_ROOT" || exit 1

echo "🔧 Reindexing entire RAG system..."
python3 rag_reindex_all.py --collection-key 3RGXZTAY
EOF

#############################
# 2. Rebuild Zotero Collection
#############################
cat > "$SCRIPT_DIR/rag_rebuild_zotero.sh" << 'EOF'
#!/bin/bash
# @raycast.title Rebuild Zotero Collection
# @raycast.author Amprenta
# @raycast.description Rebuild RAG data for the configured Zotero collection.
# @raycast.icon 📚
# @raycast.mode fullOutput
# @raycast.package Amprenta RAG
# @raycast.refreshTime 1h

RAG_ROOT="$HOME/Documents/Notion RAG"
cd "$RAG_ROOT" || exit 1

echo "📚 Rebuilding Zotero collection..."
python3 rebuild_collection_universe.py --collection-key 3RGXZTAY
EOF

#############################
# 3. Ingest Emails & Notes
#############################
cat > "$SCRIPT_DIR/rag_ingest_emails.sh" << 'EOF'
#!/bin/bash
# @raycast.title Ingest Emails & Notes
# @raycast.author Amprenta
# @raycast.description Ingest new items from the Email & Notes Inbox into the RAG Engine.
# @raycast.icon 📥
# @raycast.mode fullOutput
# @raycast.package Amprenta RAG
# @raycast.refreshTime 1h

RAG_ROOT="$HOME/Documents/Notion RAG"
cd "$RAG_ROOT" || exit 1

echo "📥 Ingesting emails & notes..."
python3 email_ingestion.py
EOF

#############################
# 4. Verify RAG Metadata
#############################
cat > "$SCRIPT_DIR/rag_verify.sh" << 'EOF'
#!/bin/bash
# @raycast.title Verify RAG Metadata
# @raycast.author Amprenta
# @raycast.description Verify consistency between Literature/Email DBs and RAG Engine chunks.
# @raycast.icon 📈
# @raycast.mode fullOutput
# @raycast.package Amprenta RAG
# @raycast.refreshTime 1h

RAG_ROOT="$HOME/Documents/Notion RAG"
cd "$RAG_ROOT" || exit 1

echo "📈 Verifying RAG metadata..."
python3 rag_verify_metadata.py
EOF

#############################
# 5. Query RAG Engine (with Raycast prompt)
#############################
cat > "$SCRIPT_DIR/rag_query.sh" << 'EOF'
#!/bin/bash
# @raycast.title Query RAG Engine
# @raycast.author Amprenta
# @raycast.description Ask a question against the RAG index (papers + emails + notes).
# @raycast.icon 🧠
# @raycast.mode fullOutput
# @raycast.package Amprenta RAG
# @raycast.argument1 label="Query" placeholder="ceramide metabolism"

RAG_ROOT="$HOME/Documents/Notion RAG"
cd "$RAG_ROOT" || exit 1

QUERY="$1"

if [ -z "$QUERY" ]; then
  echo "❌ No query provided."
  exit 1
fi

echo "🔎 Querying RAG Engine: $QUERY"
python3 rag_query.py --query "$QUERY" --top-k 10
EOF

echo "Making scripts executable..."
chmod +x "$SCRIPT_DIR"/*.sh

echo "✅ Raycast script commands installed in: $SCRIPT_DIR"
echo "➡ Now open Raycast Settings → Extensions → Script Commands → Add Script Directory…, and select:"
echo "   $SCRIPT_DIR"
